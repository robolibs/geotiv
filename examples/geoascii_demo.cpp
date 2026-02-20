#include <datapod/datapod.hpp>
namespace dp = datapod;
#include "rastkit/rastkit.hpp"
#include <iostream>

int main() {
    size_t rows = 100, cols = 100;
    double cellSize = 1.0;
    dp::Geo datum{48.0, 11.0, 500.0};
    dp::Pose shift{dp::Point{0, 0, 0}, dp::Quaternion::from_euler(0, 0, 0)};

    auto grid = dp::make_grid<float>(rows, cols, cellSize, true, shift, 0.0f);

    rastkit::RasterCollection rc;
    rc.datum = datum;
    rc.shift = shift;
    rc.resolution = cellSize;

    rastkit::Layer layer;
    layer.grid = std::move(grid);
    layer.width = static_cast<uint32_t>(cols);
    layer.height = static_cast<uint32_t>(rows);
    layer.samplesPerPixel = 1;
    layer.planarConfig = 1;
    layer.datum = datum;
    layer.shift = shift;
    layer.resolution = cellSize;
    layer.geoAsciiParams = "WGS 84 / World Geodetic System 1984";

    rc.layers.push_back(std::move(layer));

    rastkit::WriteRasterCollection(rc, "geoascii_demo.tif");

    std::cout << "Created geoascii_demo.tif with GeoAsciiParams" << std::endl;
    std::cout << "Run: gdalinfo geoascii_demo.tif" << std::endl;

    return 0;
}
