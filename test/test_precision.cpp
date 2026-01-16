#include <datapod/datapod.hpp>
namespace dp = datapod;
#include "geotiv/geotiv.hpp"
#include <cmath>
#include <doctest/doctest.h>
#include <filesystem>
#include <fstream>
#include <iomanip>

TEST_CASE("Coordinate precision at different latitudes") {
    // Test precision of round-trip conversions (save -> load) at various latitudes
    // Goal: sub-meter precision everywhere, sub-centimeter at mid-latitudes

    auto test_precision_at_location = [](const char *location, double lat, double lon, double expected_error_m) {
        SUBCASE(location) {
            const size_t rows = 10, cols = 10;
            const double cellSize = 1.0; // 1 meter resolution

            dp::Geo datum{lat, lon, 100.0};
            auto rotation = dp::Quaternion::from_euler(0.0, 0.0, 0.0);
            dp::Pose shift{dp::Point{50.0, 50.0, 10.0}, rotation};

            // Create grid
            auto grid = dp::make_grid<uint8_t>(rows, cols, cellSize, true, shift, uint8_t{128});

            // Create RasterCollection
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

            rc.layers.push_back(std::move(layer));

            // Write to file
            std::string testFile = std::string("test_precision_") + location + ".tif";
            geotiv::WriteRasterCollection(rc, testFile);

            // Read back
            auto rc2 = geotiv::ReadRasterCollection(testFile);
            REQUIRE(rc2.layers.size() == 1);

            const auto &grid2 = rc2.layers[0].gridAs<uint8_t>();

            // Test multiple points across the grid
            double max_error = 0.0;
            double sum_error = 0.0;
            int num_points = 0;

            for (size_t r = 0; r < rows; r += 3) {
                for (size_t c = 0; c < cols; c += 3) {
                    auto orig = rc.layers[0].gridAs<uint8_t>().get_point(r, c);
                    auto recon = grid2.get_point(r, c);

                    double dx = recon.x - orig.x;
                    double dy = recon.y - orig.y;
                    double dz = recon.z - orig.z;
                    double error = std::sqrt(dx * dx + dy * dy + dz * dz);

                    max_error = std::max(max_error, error);
                    sum_error += error;
                    num_points++;
                }
            }

            double avg_error = sum_error / num_points;

            // Log results
            MESSAGE("Location: ", location);
            MESSAGE("Latitude: ", lat, "°");
            MESSAGE("Max error: ", max_error, " m");
            MESSAGE("Avg error: ", avg_error, " m");

            // Check precision
            CHECK(max_error < expected_error_m);
            CHECK(avg_error < expected_error_m / 2.0);

            // Verify resolution is preserved
            double resolution_error = std::abs(rc2.layers[0].resolution - cellSize);
            CHECK(resolution_error < 0.001); // 1mm resolution error

            // Clean up
            std::filesystem::remove(testFile);
        }
    };

    // Test at different latitudes
    test_precision_at_location("Equator", 0.0, 0.0, 0.01);        // 1cm at equator
    test_precision_at_location("Mid_Latitude", 45.0, 0.0, 0.01);  // 1cm at 45°
    test_precision_at_location("High_Latitude", 70.0, 0.0, 0.01); // 1cm at 70°
    test_precision_at_location("Near_Pole", 85.0, 0.0, 0.05);     // 5cm near pole (acceptable)
}

TEST_CASE("Precision with large ENU offsets") {
    // Test that precision is maintained even with large offsets from datum
    SUBCASE("1km offset") {
        const size_t rows = 5, cols = 5;
        const double cellSize = 0.5; // 50cm resolution

        dp::Geo datum{45.0, 9.0, 100.0};
        auto rotation = dp::Quaternion::from_euler(0.0, 0.0, 0.0);
        dp::Pose shift{dp::Point{1000.0, 1000.0, 50.0}, rotation}; // 1km offset

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

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_precision_large_offset.tif";
        geotiv::WriteRasterCollection(rc, testFile);

        auto rc2 = geotiv::ReadRasterCollection(testFile);
        const auto &grid2 = rc2.layers[0].gridAs<uint8_t>();

        // Check corner points
        auto orig_00 = rc.layers[0].gridAs<uint8_t>().get_point(0, 0);
        auto recon_00 = grid2.get_point(0, 0);

        double dx = recon_00.x - orig_00.x;
        double dy = recon_00.y - orig_00.y;
        double dz = recon_00.z - orig_00.z;
        double error = std::sqrt(dx * dx + dy * dy + dz * dz);

        MESSAGE("1km offset error: ", error, " m");
        CHECK(error < 0.01); // Still sub-centimeter

        std::filesystem::remove(testFile);
    }

    SUBCASE("10km offset") {
        const size_t rows = 5, cols = 5;
        const double cellSize = 1.0;

        dp::Geo datum{45.0, 9.0, 100.0};
        auto rotation = dp::Quaternion::from_euler(0.0, 0.0, 0.0);
        dp::Pose shift{dp::Point{10000.0, 10000.0, 50.0}, rotation}; // 10km offset

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

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_precision_10km_offset.tif";
        geotiv::WriteRasterCollection(rc, testFile);

        auto rc2 = geotiv::ReadRasterCollection(testFile);
        const auto &grid2 = rc2.layers[0].gridAs<uint8_t>();

        auto orig_00 = rc.layers[0].gridAs<uint8_t>().get_point(0, 0);
        auto recon_00 = grid2.get_point(0, 0);

        double dx = recon_00.x - orig_00.x;
        double dy = recon_00.y - orig_00.y;
        double dz = recon_00.z - orig_00.z;
        double error = std::sqrt(dx * dx + dy * dy + dz * dz);

        MESSAGE("10km offset error: ", error, " m");
        CHECK(error < 0.05); // 5cm acceptable at 10km

        std::filesystem::remove(testFile);
    }
}

TEST_CASE("Precision with fine resolutions") {
    // Test sub-meter resolutions (e.g., 10cm, 1cm)
    SUBCASE("10cm resolution") {
        const size_t rows = 20, cols = 20;
        const double cellSize = 0.1; // 10cm

        dp::Geo datum{48.0, 11.0, 500.0};
        auto rotation = dp::Quaternion::from_euler(0.0, 0.0, 0.0);
        dp::Pose shift{dp::Point{10.0, 10.0, 5.0}, rotation};

        auto grid = dp::make_grid<uint16_t>(rows, cols, cellSize, true, shift, uint16_t{1000});

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

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_precision_10cm.tif";
        geotiv::WriteRasterCollection(rc, testFile);

        auto rc2 = geotiv::ReadRasterCollection(testFile);

        // Check resolution preservation
        double resolution_error = std::abs(rc2.layers[0].resolution - cellSize);
        MESSAGE("10cm resolution error: ", resolution_error, " m");
        CHECK(resolution_error < 0.0001); // 0.1mm resolution error

        const auto &grid2 = rc2.layers[0].gridAs<uint16_t>();

        // Check spatial precision
        auto orig = rc.layers[0].gridAs<uint16_t>().get_point(10, 10);
        auto recon = grid2.get_point(10, 10);

        double dx = recon.x - orig.x;
        double dy = recon.y - orig.y;
        double error = std::sqrt(dx * dx + dy * dy);

        MESSAGE("10cm grid spatial error: ", error, " m");
        CHECK(error < 0.001); // 1mm spatial error

        std::filesystem::remove(testFile);
    }
}

TEST_CASE("Precision with rotation") {
    // Test that rotation doesn't degrade precision
    SUBCASE("45 degree rotation") {
        const size_t rows = 10, cols = 10;
        const double cellSize = 1.0;

        dp::Geo datum{45.0, 9.0, 100.0};
        auto rotation = dp::Quaternion::from_euler(0.0, 0.0, M_PI / 4.0); // 45° yaw
        dp::Pose shift{dp::Point{100.0, 100.0, 50.0}, rotation};

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

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_precision_rotation.tif";
        geotiv::WriteRasterCollection(rc, testFile);

        auto rc2 = geotiv::ReadRasterCollection(testFile);
        const auto &grid2 = rc2.layers[0].gridAs<uint8_t>();

        // Test multiple points
        double max_error = 0.0;
        for (size_t r = 0; r < rows; r += 3) {
            for (size_t c = 0; c < cols; c += 3) {
                auto orig = rc.layers[0].gridAs<uint8_t>().get_point(r, c);
                auto recon = grid2.get_point(r, c);

                double dx = recon.x - orig.x;
                double dy = recon.y - orig.y;
                double dz = recon.z - orig.z;
                double error = std::sqrt(dx * dx + dy * dy + dz * dz);

                max_error = std::max(max_error, error);
            }
        }

        MESSAGE("Rotated grid max error: ", max_error, " m");
        CHECK(max_error < 0.01); // Sub-centimeter even with rotation

        std::filesystem::remove(testFile);
    }
}

TEST_CASE("Resolution conversion precision") {
    // Test that meter <-> degree conversions are precise
    SUBCASE("Various resolutions at mid-latitude") {
        dp::Geo datum{45.0, 9.0, 100.0};
        auto rotation = dp::Quaternion::from_euler(0.0, 0.0, 0.0);
        dp::Pose shift{dp::Point{0.0, 0.0, 0.0}, rotation};

        std::vector<double> resolutions = {0.01, 0.1, 0.5, 1.0, 2.5, 5.0, 10.0, 100.0};

        for (double res : resolutions) {
            auto grid = dp::make_grid<uint8_t>(5, 5, res, true, shift, uint8_t{0});

            geotiv::RasterCollection rc;
            rc.datum = datum;
            rc.shift = shift;
            rc.resolution = res;

            geotiv::Layer layer;
            layer.grid = std::move(grid);
            layer.width = 5;
            layer.height = 5;
            layer.samplesPerPixel = 1;
            layer.planarConfig = 1;
            layer.datum = datum;
            layer.shift = shift;
            layer.resolution = res;

            rc.layers.push_back(std::move(layer));

            std::string testFile = "test_resolution_" + std::to_string(res) + ".tif";
            geotiv::WriteRasterCollection(rc, testFile);

            auto rc2 = geotiv::ReadRasterCollection(testFile);

            double resolution_error = std::abs(rc2.layers[0].resolution - res);
            double relative_error = resolution_error / res;

            MESSAGE("Resolution: ", res, "m, Error: ", resolution_error, "m, Relative: ", relative_error * 100, "%");

            // Relative error should be < 0.01% for all resolutions
            CHECK(relative_error < 0.0001);

            std::filesystem::remove(testFile);
        }
    }
}
