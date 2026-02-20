#pragma once

#include "rastkit/parser.hpp"
#include "rastkit/types.hpp"
#include "rastkit/writter.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace rastkit {

    struct GridLayer {
        dp::Grid<uint8_t> grid;
        std::string name;
        std::string type;
        std::unordered_map<std::string, std::string> properties;
        std::map<uint16_t, std::vector<uint32_t>> customTags;

        inline GridLayer(const dp::Grid<uint8_t> &g, const std::string &layer_name, const std::string &layer_type = "",
                         const std::unordered_map<std::string, std::string> &props = {})
            : grid(g), name(layer_name), type(layer_type), properties(props) {}

        /// Constructor from GridVariant - extracts uint8_t grid or converts
        /// Note: RGBA grids are converted to grayscale using luminance formula
        inline GridLayer(const GridVariant &gv, const std::string &layer_name, const std::string &layer_type = "",
                         const std::unordered_map<std::string, std::string> &props = {})
            : name(layer_name), type(layer_type), properties(props) {
            // If it's already uint8_t, just copy it
            if (auto *g8 = get_grid_if<uint8_t>(gv)) {
                grid = *g8;
            } else {
                // Convert from other types to uint8_t
                auto [rows, cols] = get_grid_dimensions(gv);
                auto resolution = get_grid_resolution(gv);
                auto pose = get_grid_pose(gv);
                grid = dp::make_grid<uint8_t>(rows, cols, resolution, true, pose, uint8_t{0});

                std::visit(
                    [this](const auto &g) {
                        using GridType = std::decay_t<decltype(g)>;
                        using T = grid_element_type_t<GridType>;

                        for (size_t r = 0; r < g.rows; ++r) {
                            for (size_t c = 0; c < g.cols; ++c) {
                                if constexpr (std::is_same_v<T, RGBA>) {
                                    // Convert RGBA to grayscale using luminance formula
                                    const RGBA &rgba = g(r, c);
                                    double lum = 0.299 * rgba.r + 0.587 * rgba.g + 0.114 * rgba.b;
                                    grid(r, c) = static_cast<uint8_t>(std::clamp(lum, 0.0, 255.0));
                                } else {
                                    // Clamp scalar types to uint8_t range
                                    auto val = static_cast<double>(g(r, c));
                                    grid(r, c) = static_cast<uint8_t>(std::clamp(val, 0.0, 255.0));
                                }
                            }
                        }
                    },
                    gv);
            }
        }

        inline void setGlobalProperty(const std::string &key, const std::string &value) {
            std::hash<std::string> hasher;
            uint16_t tag = GLOBAL_PROPERTIES_BASE_TAG + (hasher(key) % 1000);
            customTags[tag] = stringToAsciiTag(key + "=" + value);
        }

        inline std::unordered_map<std::string, std::string> getGlobalProperties() const {
            std::unordered_map<std::string, std::string> props;
            for (const auto &[tag, data] : customTags) {
                if (tag >= GLOBAL_PROPERTIES_BASE_TAG && tag < GLOBAL_PROPERTIES_BASE_TAG + 1000) {
                    std::string keyValue = asciiTagToString(data);
                    size_t eq_pos = keyValue.find('=');
                    if (eq_pos != std::string::npos) {
                        std::string key = keyValue.substr(0, eq_pos);
                        std::string value = keyValue.substr(eq_pos + 1);
                        props[key] = value;
                    }
                }
            }
            return props;
        }
    };

    class Raster {
      private:
        std::vector<GridLayer> grid_layers_;
        dp::Geo datum_;
        dp::Pose shift_;
        double resolution_;

      public:
        inline Raster(const dp::Geo &datum = dp::Geo{0.001, 0.001, 1.0}, const dp::Pose &shift = dp::Pose{},
                      double resolution = 1.0)
            : datum_(datum), shift_(shift), resolution_(resolution) {}

        inline static Raster fromFile(const std::filesystem::path &path) {
            auto rc = ReadRasterCollection(path);

            if (rc.layers.empty()) {
                throw std::runtime_error("Raster::fromFile: No layers found in file");
            }

            Raster raster(rc.datum, rc.shift, rc.resolution);

            for (const auto &layer : rc.layers) {
                std::string layerName = "layer_" + std::to_string(layer.ifdOffset);
                std::string layerType = "unknown";
                std::unordered_map<std::string, std::string> props;

                if (!layer.imageDescription.empty()) {
                    std::istringstream ss(layer.imageDescription);
                    std::string token;
                    while (ss >> token) {
                        if (token == "NAME" && ss >> token) {
                            layerName = token;
                        } else if (token == "TYPE" && ss >> token) {
                            layerType = token;
                        }
                    }
                    props["description"] = layer.imageDescription;
                }

                props["width"] = std::to_string(layer.width);
                props["height"] = std::to_string(layer.height);
                props["resolution"] = std::to_string(layer.resolution);
                props["samples_per_pixel"] = std::to_string(layer.samplesPerPixel);

                GridLayer gridLayer(layer.grid, layerName, layerType, props);
                gridLayer.customTags = layer.customTags;

                raster.grid_layers_.push_back(std::move(gridLayer));
            }

            return raster;
        }

        inline void toFile(const std::filesystem::path &path) const {
            RasterCollection rc;
            rc.datum = datum_;
            rc.shift = shift_;
            rc.resolution = resolution_;

            for (const auto &gridLayer : grid_layers_) {
                Layer layer;
                layer.grid = gridLayer.grid;
                layer.width = static_cast<uint32_t>(gridLayer.grid.cols);
                layer.height = static_cast<uint32_t>(gridLayer.grid.rows);
                layer.resolution = gridLayer.grid.resolution;
                layer.datum = datum_;
                layer.shift = shift_;
                layer.samplesPerPixel = 1;
                layer.planarConfig = 1;
                layer.imageDescription = "";
                layer.customTags = gridLayer.customTags;

                rc.layers.push_back(std::move(layer));
            }

            WriteRasterCollection(rc, path);
        }

        inline size_t gridCount() const { return grid_layers_.size(); }
        inline bool hasGrids() const { return !grid_layers_.empty(); }
        inline void clearGrids() { grid_layers_.clear(); }

        inline const GridLayer &getGrid(size_t index) const {
            if (index >= grid_layers_.size()) {
                throw std::out_of_range("Grid index out of range");
            }
            return grid_layers_[index];
        }

        inline GridLayer &getGrid(size_t index) {
            if (index >= grid_layers_.size()) {
                throw std::out_of_range("Grid index out of range");
            }
            return grid_layers_[index];
        }

        inline const GridLayer &getGrid(const std::string &name) const {
            auto it = std::find_if(grid_layers_.begin(), grid_layers_.end(),
                                   [&name](const GridLayer &layer) { return layer.name == name; });
            if (it == grid_layers_.end()) {
                throw std::runtime_error("Grid with name '" + name + "' not found");
            }
            return *it;
        }

        inline GridLayer &getGrid(const std::string &name) {
            auto it = std::find_if(grid_layers_.begin(), grid_layers_.end(),
                                   [&name](const GridLayer &layer) { return layer.name == name; });
            if (it == grid_layers_.end()) {
                throw std::runtime_error("Grid with name '" + name + "' not found");
            }
            return *it;
        }

        inline void addGrid(uint32_t width, uint32_t height, const std::string &name, const std::string &type = "",
                            const std::unordered_map<std::string, std::string> &properties = {}) {
            auto grid = dp::make_grid<uint8_t>(height, width, resolution_, true, shift_, uint8_t{0});
            auto props = properties;
            if (!type.empty()) {
                props["type"] = type;
            }
            grid_layers_.emplace_back(grid, name, type, props);

            if (!grid_layers_.empty() && grid_layers_.size() > 1) {
                auto globalProps = grid_layers_[0].getGlobalProperties();
                for (const auto &[key, value] : globalProps) {
                    grid_layers_.back().setGlobalProperty(key, value);
                }
            }
        }

        inline void removeGrid(size_t index) {
            if (index < grid_layers_.size()) {
                grid_layers_.erase(grid_layers_.begin() + index);
            }
        }

        inline void addTerrainGrid(uint32_t width, uint32_t height, const std::string &name = "terrain") {
            addGrid(width, height, name, "terrain");
        }

        inline void addOcclusionGrid(uint32_t width, uint32_t height, const std::string &name = "occlusion") {
            addGrid(width, height, name, "occlusion");
        }

        inline void addElevationGrid(uint32_t width, uint32_t height, const std::string &name = "elevation") {
            addGrid(width, height, name, "elevation");
        }

        inline std::vector<GridLayer> getGridsByType(const std::string &type) const {
            std::vector<GridLayer> result;
            for (const auto &layer : grid_layers_) {
                if (layer.type == type) {
                    result.push_back(layer);
                }
            }
            return result;
        }

        inline std::vector<GridLayer> filterByProperty(const std::string &key, const std::string &value) const {
            std::vector<GridLayer> result;
            for (const auto &layer : grid_layers_) {
                auto it = layer.properties.find(key);
                if (it != layer.properties.end() && it->second == value) {
                    result.push_back(layer);
                }
            }
            return result;
        }

        inline std::vector<std::string> getGridNames() const {
            std::vector<std::string> names;
            for (const auto &layer : grid_layers_) {
                names.push_back(layer.name);
            }
            return names;
        }

        inline const dp::Geo &getDatum() const { return datum_; }
        inline void setDatum(const dp::Geo &datum) { datum_ = datum; }

        inline const dp::Pose &getShift() const { return shift_; }
        inline void setShift(const dp::Pose &shift) { shift_ = shift; }

        inline double getResolution() const { return resolution_; }
        inline void setResolution(double resolution) { resolution_ = resolution; }

        inline void setGlobalProperty(const std::string &key, const std::string &value) {
            for (auto &layer : grid_layers_) {
                layer.setGlobalProperty(key, value);
            }
        }

        inline std::string getGlobalProperty(const std::string &key, const std::string &default_value = "") const {
            if (!grid_layers_.empty()) {
                auto props = grid_layers_[0].getGlobalProperties();
                auto it = props.find(key);
                return (it != props.end()) ? it->second : default_value;
            }
            return default_value;
        }

        inline std::unordered_map<std::string, std::string> getGlobalProperties() const {
            return grid_layers_.empty() ? std::unordered_map<std::string, std::string>{}
                                        : grid_layers_[0].getGlobalProperties();
        }

        inline void removeGlobalProperty(const std::string &key) {
            std::hash<std::string> hasher;
            uint16_t tag = GLOBAL_PROPERTIES_BASE_TAG + (hasher(key) % 1000);
            for (auto &layer : grid_layers_) {
                layer.customTags.erase(tag);
            }
        }

        auto begin() { return grid_layers_.begin(); }
        auto end() { return grid_layers_.end(); }
        auto begin() const { return grid_layers_.begin(); }
        auto end() const { return grid_layers_.end(); }
        auto cbegin() const { return grid_layers_.cbegin(); }
        auto cend() const { return grid_layers_.cend(); }
    };

    template <typename T = uint8_t> dp::Layer<T> ReadLayerCollection(const std::filesystem::path &path) {
        auto rc = ReadRasterCollection(path);

        if (rc.layers.empty()) {
            throw std::runtime_error("No layers found in file");
        }

        double resolution = rc.resolution;
        dp::Pose shift = rc.shift;

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

        // Create layer using factory function
        auto layer3d = dp::make_layer<T>(rows, cols, layerCount, resolution, layerHeight, true, shift, T{});

        for (size_t layerIdx = 0; layerIdx < layerCount; ++layerIdx) {
            const auto &ifdLayer = rc.layers[layerIdx];

            // Use std::visit to handle the grid variant
            std::visit(
                [&](const auto &grid2d) {
                    using GridType = std::decay_t<decltype(grid2d)>;
                    using SrcT = grid_element_type_t<GridType>;

                    for (size_t r = 0; r < rows; ++r) {
                        for (size_t c = 0; c < cols; ++c) {
                            if constexpr (std::is_same_v<SrcT, RGBA>) {
                                // Convert RGBA to scalar using luminance formula
                                const RGBA &rgba = grid2d(r, c);
                                double lum = 0.299 * rgba.r + 0.587 * rgba.g + 0.114 * rgba.b;
                                layer3d(r, c, layerIdx) = static_cast<T>(lum);
                            } else {
                                layer3d(r, c, layerIdx) = static_cast<T>(grid2d(r, c));
                            }
                        }
                    }
                },
                ifdLayer.grid);
        }

        return layer3d;
    }

} // namespace rastkit
