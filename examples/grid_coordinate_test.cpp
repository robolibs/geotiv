#include "rastkit/rastkit.hpp"
#include <datapod/datapod.hpp>
#include <iostream>

namespace dp = datapod;

int main() {
    std::cout << "=== Grid Coordinate Test ===" << std::endl;

    // Create a Grid with the same parameters as the Layer demo
    const size_t rows = 6, cols = 8;
    const double cellSize = 2.5;

    dp::Geo datum{47.6062, -122.3321, 56.0}; // Seattle coordinates
    auto rotation = dp::Quaternion::from_euler(0.0, 0.0, 0.1);
    dp::Pose shift{dp::Point{100.0, 200.0, 50.0}, rotation}; // Same shift as Layer demo

    std::cout << "Creating 2D Grid: " << rows << "x" << cols << std::endl;
    std::cout << "Cell size: " << cellSize << "m" << std::endl;
    std::cout << "Position: (" << shift.point.x << ", " << shift.point.y << ", " << shift.point.z << ")" << std::endl;

    // Create the 2D grid
    auto originalGrid = dp::make_grid<uint8_t>(rows, cols, cellSize, true, shift, uint8_t{0});

    // Fill with test data
    for (size_t r = 0; r < rows; ++r) {
        for (size_t c = 0; c < cols; ++c) {
            originalGrid(r, c) = static_cast<uint8_t>(r * 10 + c);
        }
    }

    // Print original coordinates
    std::cout << "\nOriginal Grid coordinates:" << std::endl;
    auto origPoint00 = originalGrid.get_point(0, 0);
    auto origPoint01 = originalGrid.get_point(0, 1);
    auto origPoint10 = originalGrid.get_point(1, 0);
    std::cout << "Point (0,0): (" << origPoint00.x << ", " << origPoint00.y << ", " << origPoint00.z << ")"
              << std::endl;
    std::cout << "Point (0,1): (" << origPoint01.x << ", " << origPoint01.y << ", " << origPoint01.z << ")"
              << std::endl;
    std::cout << "Point (1,0): (" << origPoint10.x << ", " << origPoint10.y << ", " << origPoint10.z << ")"
              << std::endl;

    // Save to GeoTIFF using the existing grid functionality
    const std::string filename = "grid_coordinate_test.tif";
    std::cout << "\nSaving Grid to GeoTIFF: " << filename << std::endl;

    try {
        // Create RasterCollection with single layer
        rastkit::RasterCollection rc;
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = cellSize;

        rastkit::Layer layer;
        layer.grid = originalGrid;
        layer.width = static_cast<uint32_t>(cols);
        layer.height = static_cast<uint32_t>(rows);
        layer.samplesPerPixel = 1;
        layer.planarConfig = 1;
        layer.datum = datum;
        layer.shift = shift;
        layer.resolution = cellSize;
        // Don't set custom imageDescription - let it auto-generate with geospatial info

        rc.layers.push_back(std::move(layer));

        rastkit::WriteRasterCollection(rc, filename);
        std::cout << "Successfully saved Grid to " << filename << std::endl;

    } catch (const std::exception &e) {
        std::cerr << "Error saving: " << e.what() << std::endl;
        return 1;
    }

    // Read back and check coordinates
    std::cout << "\nReading back from GeoTIFF..." << std::endl;

    try {
        auto rc = rastkit::ReadRasterCollection(filename);
        std::cout << "Successfully read from " << filename << std::endl;

        if (rc.layers.empty()) {
            std::cerr << "No layers found!" << std::endl;
            return 1;
        }

        // Get the grid as uint8_t (we know it's uint8_t since we wrote it that way)
        const auto &reconstructedGrid = rc.layers[0].gridAs<uint8_t>();

        std::cout << "\nReconstructed Grid info:" << std::endl;
        std::cout << "Dimensions: " << reconstructedGrid.rows << "x" << reconstructedGrid.cols << std::endl;
        std::cout << "Cell size: " << reconstructedGrid.resolution << "m" << std::endl;
        std::cout << "RC shift: (" << rc.shift.point.x << ", " << rc.shift.point.y << ", " << rc.shift.point.z << ")"
                  << std::endl;
        std::cout << "RC resolution: " << rc.resolution << std::endl;

        // Check reconstructed coordinates
        std::cout << "\nReconstructed Grid coordinates:" << std::endl;
        auto reconPoint00 = reconstructedGrid.get_point(0, 0);
        auto reconPoint01 = reconstructedGrid.get_point(0, 1);
        auto reconPoint10 = reconstructedGrid.get_point(1, 0);
        std::cout << "Point (0,0): (" << reconPoint00.x << ", " << reconPoint00.y << ", " << reconPoint00.z << ")"
                  << std::endl;
        std::cout << "Point (0,1): (" << reconPoint01.x << ", " << reconPoint01.y << ", " << reconPoint01.z << ")"
                  << std::endl;
        std::cout << "Point (1,0): (" << reconPoint10.x << ", " << reconPoint10.y << ", " << reconPoint10.z << ")"
                  << std::endl;

        // Compare coordinates
        std::cout << "\nCoordinate comparison:" << std::endl;

        auto checkPoint = [](const char *name, dp::Point orig, dp::Point recon) {
            bool matches = (std::abs(orig.x - recon.x) < 0.001 && std::abs(orig.y - recon.y) < 0.001 &&
                            std::abs(orig.z - recon.z) < 0.001);
            std::cout << name << ": " << (matches ? "MATCH" : "MISMATCH") << std::endl;
            if (!matches) {
                std::cout << "  Original: (" << orig.x << ", " << orig.y << ", " << orig.z << ")" << std::endl;
                std::cout << "  Reconstructed: (" << recon.x << ", " << recon.y << ", " << recon.z << ")" << std::endl;
                std::cout << "  Difference: (" << (recon.x - orig.x) << ", " << (recon.y - orig.y) << ", "
                          << (recon.z - orig.z) << ")" << std::endl;
            }
        };

        checkPoint("Point (0,0)", origPoint00, reconPoint00);
        checkPoint("Point (0,1)", origPoint01, reconPoint01);
        checkPoint("Point (1,0)", origPoint10, reconPoint10);

        // Check data integrity too
        std::cout << "\nData integrity check:" << std::endl;
        bool dataMatches = true;
        for (size_t r = 0; r < rows; ++r) {
            for (size_t c = 0; c < cols; ++c) {
                if (originalGrid(r, c) != reconstructedGrid(r, c)) {
                    dataMatches = false;
                    std::cout << "Data mismatch at (" << r << "," << c << "): " << (int)originalGrid(r, c) << " vs "
                              << (int)reconstructedGrid(r, c) << std::endl;
                    break;
                }
            }
            if (!dataMatches)
                break;
        }

        if (dataMatches) {
            std::cout << "All data matches!" << std::endl;
        } else {
            std::cout << "Data mismatch found!" << std::endl;
        }

    } catch (const std::exception &e) {
        std::cerr << "Error reading: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
