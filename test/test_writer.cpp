#include <datapod/datapod.hpp>
namespace dp = datapod;
#include "geotiv/geotiv.hpp"
#include <cstring>
#include <doctest/doctest.h>
#include <filesystem>
#include <fstream>

TEST_CASE("GeoTIFF Writer functionality") {
    SUBCASE("Create simple raster collection and convert to bytes") {
        // Create a simple 3x2 grid with checkerboard pattern
        size_t rows = 2, cols = 3;
        double cellSize = 1.0;
        dp::Geo datum{0.0, 0.0, 0.0};
        auto rotation = dp::Quaternion::from_euler(0, 0, 0);
        dp::Pose shift{dp::Point{0, 0, 0}, rotation};

        auto grid = dp::make_grid<uint8_t>(rows, cols, cellSize, true, shift, uint8_t{0});

        // Fill with simple pattern
        grid(0, 0) = 255;
        grid(0, 1) = 0;
        grid(0, 2) = 255;
        grid(1, 0) = 0;
        grid(1, 1) = 255;
        grid(1, 2) = 0;

        // Build RasterCollection
        geotiv::RasterCollection rc;
        // CRS is always WGS84
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = cellSize;

        geotiv::Layer layer;
        layer.grid = std::move(grid);
        layer.width = static_cast<uint32_t>(cols);
        layer.height = static_cast<uint32_t>(rows);
        layer.samplesPerPixel = 1;
        layer.planarConfig = 1;
        // Set per-layer metadata
        // CRS is always WGS84
        layer.datum = datum;
        layer.shift = shift;
        layer.resolution = cellSize;

        rc.layers.push_back(std::move(layer));

        // Convert to TIFF bytes
        REQUIRE_NOTHROW(geotiv::toTiffBytes(rc));

        auto bytes = geotiv::toTiffBytes(rc);
        CHECK(bytes.size() > 0);
        CHECK(bytes.size() > 8); // At least the TIFF header
    }

    SUBCASE("Empty RasterCollection should throw") {
        geotiv::RasterCollection rc;
        CHECK_THROWS_WITH(geotiv::toTiffBytes(rc), "toTiffBytes(): no layers");
    }

    SUBCASE("Write and verify file creation") {
        // Create a simple raster collection
        size_t rows = 5, cols = 5;
        double cellSize = 1.0;
        dp::Geo datum{0.0, 0.0, 0.0};
        auto rotation = dp::Quaternion::from_euler(0, 0, 0);
        dp::Pose shift{dp::Point{0, 0, 0}, rotation};

        auto grid = dp::make_grid<uint8_t>(rows, cols, cellSize, true, shift, uint8_t{0});

        // Fill with gradient
        for (size_t r = 0; r < rows; ++r) {
            for (size_t c = 0; c < cols; ++c) {
                grid(r, c) = static_cast<uint8_t>((r * cols + c) * 10);
            }
        }

        geotiv::RasterCollection rc;
        // CRS is always WGS84
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = cellSize;

        geotiv::Layer layer;
        layer.grid = std::move(grid);
        layer.width = static_cast<uint32_t>(cols);
        layer.height = static_cast<uint32_t>(rows);
        layer.samplesPerPixel = 1;
        layer.planarConfig = 1;
        // Set per-layer metadata
        // CRS is always WGS84
        layer.datum = datum;
        layer.shift = shift;
        layer.resolution = cellSize;

        rc.layers.push_back(std::move(layer));

        // Write to file
        std::string testFile = "test_output.tif";
        REQUIRE_NOTHROW(geotiv::WriteRasterCollection(rc, testFile));

        // Verify file exists and has content
        CHECK(std::filesystem::exists(testFile));
        CHECK(std::filesystem::file_size(testFile) > 0);

        // Clean up
        std::filesystem::remove(testFile);
    }

    SUBCASE("Multi-layer GeoTIFF") {
        size_t rows = 3, cols = 3;
        double cellSize = 2.0;
        dp::Geo datum{45.0, 9.0, 100.0};
        auto rotation = dp::Quaternion::from_euler(0, 0, 0);
        dp::Pose shift{dp::Point{0, 0, 0}, rotation};

        geotiv::RasterCollection rc;
        // CRS is always WGS84
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = cellSize;

        // Create two layers with different patterns
        for (int layerIdx = 0; layerIdx < 2; ++layerIdx) {
            auto grid = dp::make_grid<uint8_t>(rows, cols, cellSize, true, shift, uint8_t{0});

            for (size_t r = 0; r < rows; ++r) {
                for (size_t c = 0; c < cols; ++c) {
                    grid(r, c) = static_cast<uint8_t>((layerIdx + 1) * 50 + r * 10 + c);
                }
            }

            geotiv::Layer layer;
            layer.grid = std::move(grid);
            layer.width = static_cast<uint32_t>(cols);
            layer.height = static_cast<uint32_t>(rows);
            layer.samplesPerPixel = 1;
            layer.planarConfig = 1;
            // Set per-layer metadata - different for each layer
            // CRS is always WGS84
            layer.datum = {datum.latitude + layerIdx * 0.01, datum.longitude + layerIdx * 0.01,
                           datum.altitude + layerIdx * 10};
            layer.shift = shift;
            layer.resolution = cellSize + layerIdx * 0.1;

            rc.layers.push_back(std::move(layer));
        }

        // Should be able to convert multi-layer to bytes
        REQUIRE_NOTHROW(geotiv::toTiffBytes(rc));
        auto bytes = geotiv::toTiffBytes(rc);
        CHECK(bytes.size() > 0);

        // Write to file
        std::string testFile = "test_multilayer.tif";
        REQUIRE_NOTHROW(geotiv::WriteRasterCollection(rc, testFile));

        CHECK(std::filesystem::exists(testFile));
        CHECK(std::filesystem::file_size(testFile) > 0);

        // Clean up
        std::filesystem::remove(testFile);
    }

    SUBCASE("Tag ordering is ascending (TIFF 6.0 compliance)") {
        // Create a raster with custom tags that would be out of order if not sorted
        size_t rows = 3, cols = 3;
        double cellSize = 1.0;
        dp::Geo datum{48.0, 11.0, 500.0};
        dp::Pose shift{dp::Point{0, 0, 0}, dp::Quaternion::from_euler(0, 0, 0)};

        auto grid = dp::make_grid<uint8_t>(rows, cols, cellSize, true, shift, uint8_t{128});

        geotiv::RasterCollection rc;
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = cellSize;

        geotiv::Layer layer;
        layer.grid = std::move(grid);
        layer.width = static_cast<uint32_t>(cols);
        layer.height = static_cast<uint32_t>(rows);
        layer.samplesPerPixel = 1;
        layer.planarConfig = 1;
        layer.datum = datum;
        layer.shift = shift;
        layer.resolution = cellSize;

        // Add custom tags that would interleave with standard tags if sorted
        layer.customTags[50001] = {100};      // After GeoTIFF tags
        layer.customTags[300] = {200};        // Between standard tags (after 284, before 339)
        layer.customTags[50000] = {300, 400}; // After GeoTIFF tags (multi-value)

        rc.layers.push_back(std::move(layer));

        auto bytes = geotiv::toTiffBytes(rc);

        // Parse the IFD to verify tag ordering
        // TIFF header: bytes 0-1 = byte order, 2-3 = magic, 4-7 = IFD offset
        uint32_t ifdOffset = 0;
        std::memcpy(&ifdOffset, &bytes[4], 4);

        // Read number of entries
        uint16_t numEntries = 0;
        std::memcpy(&numEntries, &bytes[ifdOffset], 2);
        CHECK(numEntries > 0);

        // Verify tags are in ascending order
        uint16_t prevTag = 0;
        for (uint16_t e = 0; e < numEntries; ++e) {
            size_t entryOffset = ifdOffset + 2 + e * 12;
            uint16_t tag = 0;
            std::memcpy(&tag, &bytes[entryOffset], 2);

            CHECK(tag > prevTag); // Each tag must be greater than the previous
            prevTag = tag;
        }

        // Verify specific ordering: 300 should come after 284 but before 339
        // Find positions of tags 284, 300, 339
        int pos284 = -1, pos300 = -1, pos339 = -1;
        for (uint16_t e = 0; e < numEntries; ++e) {
            size_t entryOffset = ifdOffset + 2 + e * 12;
            uint16_t tag = 0;
            std::memcpy(&tag, &bytes[entryOffset], 2);

            if (tag == 284)
                pos284 = e;
            if (tag == 300)
                pos300 = e;
            if (tag == 339)
                pos339 = e;
        }

        CHECK(pos284 >= 0);
        CHECK(pos300 >= 0);
        CHECK(pos339 >= 0);
        CHECK(pos284 < pos300);
        CHECK(pos300 < pos339);
    }

    SUBCASE("Custom tags with values are correctly written and readable") {
        size_t rows = 3, cols = 3;
        double cellSize = 1.0;
        dp::Geo datum{48.0, 11.0, 500.0};
        dp::Pose shift{dp::Point{0, 0, 0}, dp::Quaternion::from_euler(0, 0, 0)};

        auto grid = dp::make_grid<uint8_t>(rows, cols, cellSize, true, shift, uint8_t{64});

        geotiv::RasterCollection rc;
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = cellSize;

        geotiv::Layer layer;
        layer.grid = std::move(grid);
        layer.width = static_cast<uint32_t>(cols);
        layer.height = static_cast<uint32_t>(rows);
        layer.samplesPerPixel = 1;
        layer.planarConfig = 1;
        layer.datum = datum;
        layer.shift = shift;
        layer.resolution = cellSize;

        // Add custom tags
        layer.customTags[50100] = {12345};
        layer.customTags[50101] = {100, 200, 300};

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_custom_tags.tif";
        geotiv::WriteRasterCollection(rc, testFile);

        // Read back and verify custom tags
        auto rc2 = geotiv::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);

        // Custom tags should be preserved
        CHECK(rc2.layers[0].customTags.count(50100) == 1);
        CHECK(rc2.layers[0].customTags[50100].size() == 1);
        CHECK(rc2.layers[0].customTags[50100][0] == 12345);

        CHECK(rc2.layers[0].customTags.count(50101) == 1);
        CHECK(rc2.layers[0].customTags[50101].size() == 3);
        CHECK(rc2.layers[0].customTags[50101][0] == 100);
        CHECK(rc2.layers[0].customTags[50101][1] == 200);
        CHECK(rc2.layers[0].customTags[50101][2] == 300);

        std::filesystem::remove(testFile);
    }
}
