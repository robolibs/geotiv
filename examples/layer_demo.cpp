#include <iostream>
#include <filesystem>
#include "concord/concord.hpp"
#include "geotiv/geotiv.hpp"
#include "geotiv/raster.hpp"

int main() {
    std::cout << "=== Layer to GeoTIFF Demo ===" << std::endl;
    
    // Create a 3D Layer with some test data
    const size_t rows = 6, cols = 8, layers = 4;
    const double cellSize = 2.5;      // 2.5 meters per cell
    const double layerHeight = 3.0;   // 3.0 meters per Z layer
    
    // Position the layer at a specific location with some rotation
    concord::Datum datum{47.6062, -122.3321, 56.0}; // Seattle coordinates
    concord::Pose shift{concord::Point{100.0, 200.0, 50.0}, concord::Euler{0.0, 0.0, 0.1}}; // 100m east, 200m north, 50m up, slight rotation
    
    std::cout << "Creating 3D Layer: " << rows << "x" << cols << "x" << layers << std::endl;
    std::cout << "Cell size: " << cellSize << "m, Layer height: " << layerHeight << "m" << std::endl;
    std::cout << "Position: (" << shift.point.x << ", " << shift.point.y << ", " << shift.point.z << ")" << std::endl;
    
    // Create the 3D layer
    concord::Layer<uint8_t> originalLayer(rows, cols, layers, cellSize, layerHeight, 
                                         true, shift, false, false);
    
    // Fill with interesting test pattern
    std::cout << "Filling with test pattern..." << std::endl;
    for (size_t l = 0; l < layers; ++l) {
        for (size_t r = 0; r < rows; ++r) {
            for (size_t c = 0; c < cols; ++c) {
                // Create a pattern that varies by layer, row, and column
                uint8_t value = static_cast<uint8_t>(
                    (l * 60) +                    // Base value per layer
                    (r * 8) +                     // Row contribution
                    (c * 2) +                     // Column contribution
                    ((r + c + l) % 3 * 20)        // Some variation
                );
                originalLayer(r, c, l) = std::min(value, uint8_t(255));
            }
        }
    }
    
    // Print some sample coordinates
    std::cout << "\nSample world coordinates:" << std::endl;
    for (size_t l = 0; l < layers; ++l) {
        auto point = originalLayer.get_point(0, 0, l);
        std::cout << "Layer " << l << " corner (0,0): (" << point.x << ", " << point.y << ", " << point.z << ")" << std::endl;
    }
    
    // Save to GeoTIFF file
    const std::string filename = "layer_demo.tif";
    std::cout << "\nSaving to GeoTIFF file: " << filename << std::endl;
    
    try {
        geotiv::WriteLayerCollection(originalLayer, filename, datum);
        std::cout << "✅ Successfully saved to " << filename << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "❌ Error saving: " << e.what() << std::endl;
        return 1;
    }
    
    // Verify file exists and show size
    if (std::filesystem::exists(filename)) {
        auto size = std::filesystem::file_size(filename);
        std::cout << "File size: " << size << " bytes" << std::endl;
    }
    
    // Read back from file
    std::cout << "\nReading back from GeoTIFF file..." << std::endl;
    
    try {
        auto reconstructedLayer = geotiv::ReadLayerCollection<uint8_t>(filename);
        std::cout << "✅ Successfully read from " << filename << std::endl;
        
        // Verify dimensions
        std::cout << "\nVerifying reconstruction:" << std::endl;
        std::cout << "Dimensions: " << reconstructedLayer.rows() << "x" 
                  << reconstructedLayer.cols() << "x" << reconstructedLayer.layers() << std::endl;
        std::cout << "Cell size: " << reconstructedLayer.inradius() << "m" << std::endl;
        std::cout << "Layer height: " << reconstructedLayer.layer_height() << "m" << std::endl;
        
        // Check coordinates match
        std::cout << "\nReconstructed coordinates:" << std::endl;
        for (size_t l = 0; l < layers; ++l) {
            auto origPoint = originalLayer.get_point(0, 0, l);
            auto reconPoint = reconstructedLayer.get_point(0, 0, l);
            std::cout << "Layer " << l << " corner (0,0): (" << reconPoint.x << ", " << reconPoint.y << ", " << reconPoint.z << ")";
            
            // Check if coordinates match
            bool matches = (std::abs(origPoint.x - reconPoint.x) < 0.001 &&
                           std::abs(origPoint.y - reconPoint.y) < 0.001 &&
                           std::abs(origPoint.z - reconPoint.z) < 0.001);
            std::cout << (matches ? " ✅" : " ❌") << std::endl;
        }
        
        // Verify data integrity
        std::cout << "\nVerifying data integrity..." << std::endl;
        bool dataMatches = true;
        size_t mismatchCount = 0;
        
        for (size_t l = 0; l < layers; ++l) {
            for (size_t r = 0; r < rows; ++r) {
                for (size_t c = 0; c < cols; ++c) {
                    if (originalLayer(r, c, l) != reconstructedLayer(r, c, l)) {
                        dataMatches = false;
                        mismatchCount++;
                        if (mismatchCount <= 3) { // Only show first few mismatches
                            std::cout << "Data mismatch at (" << r << "," << c << "," << l << "): " 
                                     << (int)originalLayer(r, c, l) << " vs " << (int)reconstructedLayer(r, c, l) << std::endl;
                        }
                    }
                }
            }
        }
        
        if (dataMatches) {
            std::cout << "✅ All data matches perfectly!" << std::endl;
        } else {
            std::cout << "❌ Found " << mismatchCount << " data mismatches" << std::endl;
        }
        
        // Show some sample data values
        std::cout << "\nSample data values (layer 0, row 0):" << std::endl;
        std::cout << "Original:      ";
        for (size_t c = 0; c < std::min(cols, size_t(8)); ++c) {
            std::cout << std::setw(4) << (int)originalLayer(0, c, 0);
        }
        std::cout << std::endl;
        
        std::cout << "Reconstructed: ";
        for (size_t c = 0; c < std::min(cols, size_t(8)); ++c) {
            std::cout << std::setw(4) << (int)reconstructedLayer(0, c, 0);
        }
        std::cout << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "❌ Error reading: " << e.what() << std::endl;
        return 1;
    }
    
    std::cout << "\n=== Demo completed successfully! ===" << std::endl;
    std::cout << "The GeoTIFF file '" << filename << "' contains the 3D layer data." << std::endl;
    std::cout << "You can inspect it with GIS tools or use it in other applications." << std::endl;
    
    return 0;
}