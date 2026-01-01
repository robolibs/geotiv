#include <datapod/datapod.hpp>
namespace dp = datapod;
#include "geotiv/geotiv.hpp"
#include <cmath>
#include <doctest/doctest.h>
#include <filesystem>

TEST_CASE("Advanced Multi-IFD with Per-Layer Tags and Metadata") {
    SUBCASE("Multi-layer with different coordinate systems and custom tags") {
        geotiv::RasterCollection rc;

        // Create 3 layers with different properties
        for (int i = 0; i < 3; ++i) {
            size_t rows = 10 + i * 5, cols = 15 + i * 5;
            double cellSize = 1.0 + i * 0.5;
            dp::Geo datum{47.0 + i * 0.1, 8.0 + i * 0.1, 100.0 + i * 50};
            double yaw_rad = i * 15.0 * M_PI / 180.0; // Different rotation per layer (in radians)
            auto rotation = dp::Quaternion::from_euler(0, 0, yaw_rad);
            dp::Pose shift{dp::Point{0, 0, 0}, rotation};

            auto grid = dp::make_grid<uint8_t>(rows, cols, cellSize, true, shift, uint8_t{0});

            // Fill with layer-specific pattern
            for (size_t r = 0; r < rows; ++r) {
                for (size_t c = 0; c < cols; ++c) {
                    grid(r, c) = static_cast<uint8_t>((i * 50 + r + c) % 256);
                }
            }

            geotiv::Layer layer;
            layer.grid = std::move(grid);
            layer.width = static_cast<uint32_t>(cols);
            layer.height = static_cast<uint32_t>(rows);
            layer.samplesPerPixel = 1;
            layer.planarConfig = 1;

            // Set different coordinate systems per layer
            // CRS is always WGS84
            layer.datum = datum;
            layer.shift = shift;
            layer.resolution = cellSize;
            // Don't set custom imageDescription - let the writer generate the proper geospatial one

            // Add custom tags specific to each layer
            layer.customTags[50000 + i] = {static_cast<uint32_t>(i * 1000)};          // Layer ID
            layer.customTags[50100] = {static_cast<uint32_t>(1735689600 + i * 3600)}; // Timestamp
            layer.customTags[50200 + i] = {42u, static_cast<uint32_t>(100 + i),
                                           static_cast<uint32_t>(200 + i)}; // Multi-value custom tag

            rc.layers.push_back(std::move(layer));
        }

        // Set collection defaults from first layer
        // CRS is always WGS84
        rc.datum = rc.layers[0].datum;
        rc.shift = rc.layers[0].shift;
        rc.resolution = rc.layers[0].resolution;

        // Write multi-IFD GeoTIFF
        std::string testFile = "multi_ifd_advanced.tif";
        REQUIRE_NOTHROW(geotiv::WriteRasterCollection(rc, testFile));
        CHECK(std::filesystem::exists(testFile));

        // Read back and verify per-IFD metadata preservation
        geotiv::RasterCollection readRc;
        REQUIRE_NOTHROW(readRc = geotiv::ReadRasterCollection(testFile));

        // Verify we have all layers
        CHECK(readRc.layers.size() == 3);

        // Verify each layer has its own metadata
        for (size_t i = 0; i < 3; ++i) {
            const auto &layer = readRc.layers[i];

            // Check dimensions
            CHECK(layer.width == (15 + i * 5));
            CHECK(layer.height == (10 + i * 5));

            // Check coordinate system
            // CRS is always WGS84
            // CRS is always WGS84

            // Check datum (with some tolerance for floating point precision)
            CHECK(layer.datum.latitude == doctest::Approx(47.0 + i * 0.1).epsilon(0.001));
            CHECK(layer.datum.longitude == doctest::Approx(8.0 + i * 0.1).epsilon(0.001));
            CHECK(layer.datum.altitude == doctest::Approx(100.0 + i * 50).epsilon(0.1));

            // Check resolution
            CHECK(layer.resolution == doctest::Approx(1.0 + i * 0.5).epsilon(0.001));

            // Check heading (in radians)
            double expected_yaw_rad = i * 15.0 * M_PI / 180.0;
            CHECK(layer.shift.rotation.to_euler().yaw == doctest::Approx(expected_yaw_rad).epsilon(0.01));

            // Check custom tags were preserved
            auto it = layer.customTags.find(50000 + i);
            CHECK(it != layer.customTags.end());
            if (it != layer.customTags.end()) {
                CHECK(it->second.size() >= 1);
                CHECK(it->second[0] == static_cast<uint32_t>(i * 1000));
            }

            // Verify grid data integrity
            auto [gridRows, gridCols] = geotiv::get_grid_dimensions(layer.grid);
            CHECK(gridRows == (10 + i * 5));
            CHECK(gridCols == (15 + i * 5));

            // Check some pixel values (grid is uint8_t)
            const auto &grid = layer.gridAs<uint8_t>();
            uint8_t expectedFirst = static_cast<uint8_t>((i * 50) % 256);
            CHECK(grid(0, 0) == expectedFirst);
        }

        // Clean up
        std::filesystem::remove(testFile);
    }

    SUBCASE("Time-series data with timestamps in custom tags") {
        geotiv::RasterCollection timeSeries;

        // Simulate a time series with 4 time points
        std::vector<uint32_t> timestamps = {1735689600, 1735693200, 1735696800, 1735700400}; // 1-hour intervals

        for (size_t t = 0; t < timestamps.size(); ++t) {
            size_t rows = 20, cols = 30;
            double cellSize = 0.5;                          // 50cm resolution
            dp::Geo surveyLocation{46.5204, 6.6234, 372.0}; // Geneva coordinates
            auto rotation = dp::Quaternion::from_euler(0, 0, 0);
            dp::Pose shift{dp::Point{0, 0, 0}, rotation};

            auto grid = dp::make_grid<uint8_t>(rows, cols, cellSize, true, shift, uint8_t{0});

            // Fill with time-dependent pattern
            for (size_t r = 0; r < rows; ++r) {
                for (size_t c = 0; c < cols; ++c) {
                    grid(r, c) = static_cast<uint8_t>((t * 40 + r + c) % 256);
                }
            }

            geotiv::Layer layer;
            layer.grid = std::move(grid);
            layer.width = static_cast<uint32_t>(cols);
            layer.height = static_cast<uint32_t>(rows);
            layer.samplesPerPixel = 1;
            layer.planarConfig = 1;
            // CRS is always WGS84
            layer.datum = surveyLocation;
            layer.shift = shift;
            layer.resolution = cellSize;
            layer.imageDescription = "Time-series data point " + std::to_string(t);

            // Store timestamp and sequence number in custom tags
            layer.customTags[50100] = {timestamps[t]};                                            // Timestamp
            layer.customTags[50101] = {static_cast<uint32_t>(t)};                                 // Sequence number
            layer.customTags[50102] = {static_cast<uint32_t>(rows), static_cast<uint32_t>(cols)}; // Dimensions

            timeSeries.layers.push_back(std::move(layer));
        }

        // Set collection metadata from first layer
        // CRS is always WGS84
        timeSeries.datum = timeSeries.layers[0].datum;
        timeSeries.shift = timeSeries.layers[0].shift;
        timeSeries.resolution = timeSeries.layers[0].resolution;

        // Write time-series GeoTIFF
        std::string timeSeriesFile = "time_series.tif";
        REQUIRE_NOTHROW(geotiv::WriteRasterCollection(timeSeries, timeSeriesFile));

        // Read back and verify temporal metadata
        geotiv::RasterCollection readTimeSeries;
        REQUIRE_NOTHROW(readTimeSeries = geotiv::ReadRasterCollection(timeSeriesFile));

        CHECK(readTimeSeries.layers.size() == 4);

        // Verify timestamps are preserved in order
        for (size_t t = 0; t < timestamps.size(); ++t) {
            const auto &layer = readTimeSeries.layers[t];

            // Check timestamp tag
            auto timestampIt = layer.customTags.find(50100);
            CHECK(timestampIt != layer.customTags.end());
            if (timestampIt != layer.customTags.end()) {
                CHECK(timestampIt->second[0] == timestamps[t]);
            }

            // Check sequence number
            auto seqIt = layer.customTags.find(50101);
            CHECK(seqIt != layer.customTags.end());
            if (seqIt != layer.customTags.end()) {
                CHECK(seqIt->second[0] == static_cast<uint32_t>(t));
            }

            // Check dimensions tag
            auto dimIt = layer.customTags.find(50102);
            CHECK(dimIt != layer.customTags.end());
            if (dimIt != layer.customTags.end() && dimIt->second.size() >= 2) {
                CHECK(dimIt->second[0] == 20); // rows
                CHECK(dimIt->second[1] == 30); // cols
            }
        }

        // Clean up
        std::filesystem::remove(timeSeriesFile);
    }
}
