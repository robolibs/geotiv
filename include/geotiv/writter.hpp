#pragma once

#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include "concord/concord.hpp"
#include "geotiv/types.hpp"

namespace geotiv {

    std::vector<uint8_t> toTiffBytes(RasterCollection const &rc);

    void WriteRasterCollection(RasterCollection const &rc, std::filesystem::path const &outPath);

    template <typename T>
    void WriteLayerCollection(const concord::Layer<T> &layer3d, const std::filesystem::path &outPath,
                              const concord::Datum &datum = concord::Datum{0.001, 0.001, 1.0}) {

        RasterCollection rc;
        rc.datum = datum;
        rc.shift = layer3d.shift();
        rc.resolution = layer3d.inradius();

        for (size_t layerIdx = 0; layerIdx < layer3d.layers(); ++layerIdx) {
            auto grid2d = layer3d.extract_grid(layerIdx);

            auto centerPoint = layer3d.get_point(0, 0, layerIdx);
            double zCoord = centerPoint.z;

            concord::Pose layerShift = layer3d.shift();
            layerShift.point.z = zCoord;
            concord::Grid<uint8_t> uint8Grid(grid2d.rows(), grid2d.cols(), grid2d.inradius(), true, layerShift);

            for (size_t r = 0; r < grid2d.rows(); ++r) {
                for (size_t c = 0; c < grid2d.cols(); ++c) {
                    uint8Grid(r, c) =
                        static_cast<uint8_t>(std::min(255.0, std::max(0.0, static_cast<double>(grid2d(r, c)))));
                }
            }

            geotiv::Layer layer;
            layer.grid = std::move(uint8Grid);
            layer.width = static_cast<uint32_t>(layer3d.cols());
            layer.height = static_cast<uint32_t>(layer3d.rows());
            layer.samplesPerPixel = 1;
            layer.planarConfig = 1;
            layer.datum = datum;
            layer.shift = layer3d.shift();
            layer.resolution = layer3d.inradius();

            layer.imageDescription =
                "CRS WGS84 DATUM " + std::to_string(datum.lat) + " " + std::to_string(datum.lon) + " " +
                std::to_string(datum.alt) + " SHIFT " + std::to_string(layer3d.shift().point.x) + " " +
                std::to_string(layer3d.shift().point.y) + " " + std::to_string(layer3d.shift().point.z) + " " +
                std::to_string(layer3d.shift().angle.yaw) + " " +
                "LayerHeight=" + std::to_string(layer3d.layer_height()) + " " +
                "Resolution=" + std::to_string(layer3d.inradius()) + " " +
                "LayerCount=" + std::to_string(layer3d.layers()) + " " + "LayerIndex=" + std::to_string(layerIdx) +
                " " + "LayerZ=" + std::to_string(zCoord);

            rc.layers.push_back(std::move(layer));
        }

        WriteRasterCollection(rc, outPath);
    }

} // namespace geotiv
