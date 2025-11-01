#pragma once

#include "geotiv.hpp"
#include <algorithm>
#include <functional>
#include <optional>
#include <stdexcept>
#include <unordered_map>

namespace geotiv {

    struct GridLayer {
        concord::Grid<uint8_t> grid;
        std::string name;
        std::string type;
        std::unordered_map<std::string, std::string> properties;
        std::map<uint16_t, std::vector<uint32_t>> customTags;

        GridLayer(const concord::Grid<uint8_t> &g, const std::string &layer_name, const std::string &layer_type = "",
                  const std::unordered_map<std::string, std::string> &props = {});

        void setGlobalProperty(const std::string &key, const std::string &value);
        std::unordered_map<std::string, std::string> getGlobalProperties() const;
    };

    class Raster {
      private:
        std::vector<GridLayer> grid_layers_;
        concord::Datum datum_;
        concord::Pose shift_;
        double resolution_;

      public:
        Raster(const concord::Datum &datum = concord::Datum{0.001, 0.001, 1.0},
               const concord::Pose &shift = concord::Pose{concord::Point{0, 0, 0}, concord::Euler{0, 0, 0}},
               double resolution = 1.0);

        static Raster fromFile(const std::filesystem::path &path);
        void toFile(const std::filesystem::path &path) const;

        size_t gridCount() const;
        bool hasGrids() const;
        void clearGrids();

        const GridLayer &getGrid(size_t index) const;
        GridLayer &getGrid(size_t index);
        const GridLayer &getGrid(const std::string &name) const;
        GridLayer &getGrid(const std::string &name);

        void addGrid(uint32_t width, uint32_t height, const std::string &name, const std::string &type = "",
                     const std::unordered_map<std::string, std::string> &properties = {});
        void removeGrid(size_t index);

        void addTerrainGrid(uint32_t width, uint32_t height, const std::string &name = "terrain");
        void addOcclusionGrid(uint32_t width, uint32_t height, const std::string &name = "occlusion");
        void addElevationGrid(uint32_t width, uint32_t height, const std::string &name = "elevation");

        std::vector<GridLayer> getGridsByType(const std::string &type) const;
        std::vector<GridLayer> filterByProperty(const std::string &key, const std::string &value) const;
        std::vector<std::string> getGridNames() const;

        const concord::Datum &getDatum() const;
        void setDatum(const concord::Datum &datum);

        const concord::Pose &getShift() const;
        void setShift(const concord::Pose &shift);

        double getResolution() const;
        void setResolution(double resolution);

        void setGlobalProperty(const std::string &key, const std::string &value);
        std::string getGlobalProperty(const std::string &key, const std::string &default_value = "") const;
        std::unordered_map<std::string, std::string> getGlobalProperties() const;
        void removeGlobalProperty(const std::string &key);

        auto begin() { return grid_layers_.begin(); }
        auto end() { return grid_layers_.end(); }
        auto begin() const { return grid_layers_.begin(); }
        auto end() const { return grid_layers_.end(); }
        auto cbegin() const { return grid_layers_.cbegin(); }
        auto cend() const { return grid_layers_.cend(); }
    };

    template <typename T = uint8_t> concord::Layer<T> ReadLayerCollection(const std::filesystem::path &path) {
        auto rc = geotiv::ReadRasterCollection(path);

        if (rc.layers.empty()) {
            throw std::runtime_error("No layers found in file");
        }

        double resolution = rc.resolution;
        concord::Pose shift = rc.shift;

        double layerHeight = 1.0;
        const auto &firstLayer = rc.layers[0];
        if (!firstLayer.imageDescription.empty()) {
            std::string desc = firstLayer.imageDescription;

            size_t pos = desc.find("LayerHeight=");
            if (pos != std::string::npos) {
                pos += 12;
                size_t endPos = desc.find(' ', pos);
                if (endPos == std::string::npos)
                    endPos = desc.length();
                std::string heightStr = desc.substr(pos, endPos - pos);
                layerHeight = std::stod(heightStr);
            }
        }

        size_t rows = firstLayer.height;
        size_t cols = firstLayer.width;
        size_t layerCount = rc.layers.size();

        concord::Pose correctShift = shift;

        concord::Layer<T> layer3d(rows, cols, layerCount, resolution, layerHeight, true, correctShift, false, false);

        for (size_t layerIdx = 0; layerIdx < layerCount; ++layerIdx) {
            const auto &ifdLayer = rc.layers[layerIdx];
            const auto &grid2d = ifdLayer.grid;

            for (size_t r = 0; r < rows; ++r) {
                for (size_t c = 0; c < cols; ++c) {
                    layer3d(r, c, layerIdx) = static_cast<T>(grid2d(r, c));
                }
            }
        }

        return layer3d;
    }

} // namespace geotiv
