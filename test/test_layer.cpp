#include <datapod/datapod.hpp>
namespace dp = datapod;
#include "rastkit/rastkit.hpp"
#include "rastkit/raster.hpp"
#include <doctest/doctest.h>

TEST_CASE("Layer to GeoTIFF conversion tests") {
    // Test parameters
    const size_t rows = 4, cols = 4, layers = 3;
    const double cellSize = 1.0;
    const double layerHeight = 2.0;
    const std::string testFile = "test_layer.tif";

    // Clean up any existing file
    std::filesystem::remove(testFile);

    SUBCASE("WriteLayerCollection and ReadLayerCollection round-trip") {
        dp::Geo datum{46.0, 8.0, 1000.0};
        auto rotation = dp::Quaternion::from_euler(0, 0, 0);
        dp::Pose shift{dp::Point{10, 20, 30}, rotation};

        // Create original 3D layer
        auto originalLayer =
            dp::make_layer<uint8_t>(rows, cols, layers, cellSize, layerHeight, true, shift, uint8_t{0});

        // Fill with test data
        for (size_t l = 0; l < layers; ++l) {
            for (size_t r = 0; r < rows; ++r) {
                for (size_t c = 0; c < cols; ++c) {
                    originalLayer(r, c, l) = static_cast<uint8_t>(l * 50 + r * 10 + c);
                }
            }
        }

        // Write to file
        REQUIRE_NOTHROW(rastkit::WriteLayerCollection(originalLayer, testFile, datum));

        // Verify file exists
        CHECK(std::filesystem::exists(testFile));

        // Read back
        dp::Layer<uint8_t> reconstructed;
        REQUIRE_NOTHROW(reconstructed = rastkit::ReadLayerCollection<uint8_t>(testFile));

        // Verify dimensions
        CHECK(reconstructed.rows == rows);
        CHECK(reconstructed.cols == cols);
        CHECK(reconstructed.layers == layers);
        CHECK(reconstructed.resolution == doctest::Approx(cellSize));
        CHECK(reconstructed.layer_height == doctest::Approx(layerHeight));

        // Verify all data
        bool dataMatches = true;
        for (size_t l = 0; l < layers; ++l) {
            for (size_t r = 0; r < rows; ++r) {
                for (size_t c = 0; c < cols; ++c) {
                    if (originalLayer(r, c, l) != reconstructed(r, c, l)) {
                        dataMatches = false;
                        INFO("Data mismatch at (" << r << "," << c << "," << l << "): " << (int)originalLayer(r, c, l)
                                                  << " vs " << (int)reconstructed(r, c, l));
                        break;
                    }
                }
                if (!dataMatches)
                    break;
            }
            if (!dataMatches)
                break;
        }
        CHECK(dataMatches == true);

        // Verify Z coordinates
        for (size_t l = 0; l < layers; ++l) {
            auto origPoint = originalLayer.get_point(0, 0, l);
            auto reconPoint = reconstructed.get_point(0, 0, l);
            CHECK(origPoint.z == doctest::Approx(reconPoint.z).epsilon(0.001));
        }

        // Clean up
        std::filesystem::remove(testFile);
    }
}
