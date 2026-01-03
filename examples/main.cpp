// main.cpp

#include <iostream>

#include "geotiv/geotiv.hpp"
#include <concord/concord.hpp>
#include <datapod/datapod.hpp>

namespace dp = datapod;

int main() {
    try {
        // --- 1) Create your Grid<uint8_t> (e.g. a 100×50 checkerboard) ---
        size_t rows = 50, cols = 100;
        double cellSize = 2.0;                               // meters per pixel
        dp::Geo datum{48.0, 11.0, 500.0};                    // lat, lon, alt
        auto rotation = dp::Quaternion::from_euler(0, 0, 0); // no rotation
        dp::Pose shift{dp::Point{0, 0, 0}, rotation};        // center grid at datum
        auto grid = dp::make_grid<uint8_t>(rows, cols, cellSize, true, shift, uint8_t{0});
        for (size_t r = 0; r < rows; ++r)
            for (size_t c = 0; c < cols; ++c)
                grid(r, c) = ((r / 5 + c / 5) % 2) ? 255 : 0;

        // --- 2) Build a RasterCollection around it ---
        geotiv::RasterCollection rc;
        // CRS is always WGS84
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = cellSize;

        geotiv::Layer L;
        L.grid = std::move(grid);
        L.width = static_cast<uint32_t>(cols);
        L.height = static_cast<uint32_t>(rows);
        L.samplesPerPixel = 1; // single‐band
        L.planarConfig = 1;    // chunky
        // Set per-layer metadata
        // CRS is always WGS84
        L.datum = datum;
        L.shift = shift;
        L.resolution = cellSize;

        rc.layers.clear();
        rc.layers.push_back(std::move(L));

        // --- 3) Write out the GeoTIFF in one call ---
        geotiv::WriteRasterCollection(rc, "output.tif");
        std::cout << "Wrote GeoTIFF " << cols << "x" << rows << " at " << cellSize << "m/px -> output.tif\n";
    } catch (std::exception &e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
