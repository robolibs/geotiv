// generate_test_pattern.cpp
// Creates a 640x640 GeoTIFF with a recognizable test pattern for image viewer compatibility testing

#include <cmath>
#include <iostream>

#include "rastkit/rastkit.hpp"
#include <concord/concord.hpp>
#include <datapod/datapod.hpp>

namespace dp = datapod;
namespace cc = concord;

int main() {
    try {
        std::cout << "Generating 640x640 test pattern GeoTIFF...\n";

        // Image dimensions
        size_t rows = 640, cols = 640;
        double cellSize = 1.0; // 1 meter per pixel

        // Use a real-world location
        dp::Geo datum{46.8182, 8.2275, 1000.0};              // lat, lon, alt
        auto rotation = dp::Quaternion::from_euler(0, 0, 0); // no rotation
        dp::Pose shift{dp::Point{0, 0, 0}, rotation};

        // Create the grid
        auto grid = dp::make_grid<uint8_t>(rows, cols, cellSize, true, shift, uint8_t{0});

        std::cout << "Creating test pattern...\n";

        // Create a recognizable test pattern
        for (size_t r = 0; r < rows; ++r) {
            for (size_t c = 0; c < cols; ++c) {
                uint8_t value = 0;

                // Create different patterns in different quadrants
                if (r < rows / 2 && c < cols / 2) {
                    // Top-left: Checkerboard pattern
                    value = ((r / 16 + c / 16) % 2) ? 255 : 64;
                } else if (r < rows / 2 && c >= cols / 2) {
                    // Top-right: Horizontal stripes
                    value = (r / 8 % 2) ? 200 : 100;
                } else if (r >= rows / 2 && c < cols / 2) {
                    // Bottom-left: Vertical stripes
                    value = (c / 8 % 2) ? 180 : 80;
                } else {
                    // Bottom-right: Concentric circles
                    double center_r = rows * 0.75;
                    double center_c = cols * 0.75;
                    double dist = std::sqrt((r - center_r) * (r - center_r) + (c - center_c) * (c - center_c));
                    value = static_cast<uint8_t>(128 + 127 * std::sin(dist / 10.0));
                }

                grid(r, c) = value;
            }
        }

        // Create RasterCollection
        rastkit::RasterCollection rc;
        // CRS is always WGS84
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = cellSize;

        // Create layer
        rastkit::Layer layer;
        layer.grid = std::move(grid);
        layer.width = static_cast<uint32_t>(cols);
        layer.height = static_cast<uint32_t>(rows);
        layer.samplesPerPixel = 1; // grayscale
        layer.planarConfig = 1;    // chunky format

        // Set per-layer geospatial metadata
        // CRS is always WGS84
        layer.datum = datum;
        layer.shift = shift;
        layer.resolution = cellSize;

        rc.layers.push_back(std::move(layer));

        // Write the GeoTIFF
        std::string filename = "test_pattern_640x640.tif";
        rastkit::WriteRasterCollection(rc, filename);

        std::cout << "Successfully created: " << filename << "\n";
        std::cout << "   Size: 640x640 pixels\n";
        std::cout << "   Type: 8-bit grayscale\n";
        std::cout << "   Format: Standard TIFF with GeoTIFF tags\n";
        std::cout << "\nPattern layout:\n";
        std::cout << "   +-------------+-------------+\n";
        std::cout << "   | Checkerboard| Horizontal  |\n";
        std::cout << "   |   pattern   |   stripes   |\n";
        std::cout << "   +-------------+-------------+\n";
        std::cout << "   |  Vertical   | Concentric  |\n";
        std::cout << "   |   stripes   |   circles   |\n";
        std::cout << "   +-------------+-------------+\n";
        std::cout << "\nThis should be viewable in any TIFF-compatible image viewer!\n";

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
