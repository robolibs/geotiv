#include "geotiv/raster.hpp"
#include <algorithm>
#include <functional>
#include <optional>
#include <stdexcept>
#include <unordered_map>

namespace geotiv {

    GridLayer::GridLayer(const concord::Grid<uint8_t> &g, const std::string &layer_name, const std::string &layer_type,
                         const std::unordered_map<std::string, std::string> &props)
        : grid(g), name(layer_name), type(layer_type), properties(props) {}

    void GridLayer::setGlobalProperty(const std::string &key, const std::string &value) {
        std::hash<std::string> hasher;
        uint16_t tag = GLOBAL_PROPERTIES_BASE_TAG + (hasher(key) % 1000);
        customTags[tag] = stringToAsciiTag(key + "=" + value);
    }

    std::unordered_map<std::string, std::string> GridLayer::getGlobalProperties() const {
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

    Raster::Raster(const concord::Datum &datum, const concord::Pose &shift, double resolution)
        : datum_(datum), shift_(shift), resolution_(resolution) {}

    Raster Raster::fromFile(const std::filesystem::path &path) {
        auto rc = geotiv::ReadRasterCollection(path);

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

    void Raster::toFile(const std::filesystem::path &path) const {
        RasterCollection rc;
        rc.datum = datum_;
        rc.shift = shift_;
        rc.resolution = resolution_;

        for (const auto &gridLayer : grid_layers_) {
            Layer layer;
            layer.grid = gridLayer.grid;
            layer.width = static_cast<uint32_t>(gridLayer.grid.cols());
            layer.height = static_cast<uint32_t>(gridLayer.grid.rows());
            layer.resolution = gridLayer.grid.inradius();
            layer.datum = datum_;
            layer.shift = shift_;
            layer.samplesPerPixel = 1;
            layer.planarConfig = 1;
            layer.imageDescription = "";
            layer.customTags = gridLayer.customTags;

            rc.layers.push_back(std::move(layer));
        }

        geotiv::WriteRasterCollection(rc, path);
    }

    size_t Raster::gridCount() const { return grid_layers_.size(); }
    bool Raster::hasGrids() const { return !grid_layers_.empty(); }
    void Raster::clearGrids() { grid_layers_.clear(); }

    const GridLayer &Raster::getGrid(size_t index) const {
        if (index >= grid_layers_.size()) {
            throw std::out_of_range("Grid index out of range");
        }
        return grid_layers_[index];
    }

    GridLayer &Raster::getGrid(size_t index) {
        if (index >= grid_layers_.size()) {
            throw std::out_of_range("Grid index out of range");
        }
        return grid_layers_[index];
    }

    const GridLayer &Raster::getGrid(const std::string &name) const {
        auto it = std::find_if(grid_layers_.begin(), grid_layers_.end(),
                               [&name](const GridLayer &layer) { return layer.name == name; });
        if (it == grid_layers_.end()) {
            throw std::runtime_error("Grid with name '" + name + "' not found");
        }
        return *it;
    }

    GridLayer &Raster::getGrid(const std::string &name) {
        auto it = std::find_if(grid_layers_.begin(), grid_layers_.end(),
                               [&name](const GridLayer &layer) { return layer.name == name; });
        if (it == grid_layers_.end()) {
            throw std::runtime_error("Grid with name '" + name + "' not found");
        }
        return *it;
    }

    void Raster::addGrid(uint32_t width, uint32_t height, const std::string &name, const std::string &type,
                         const std::unordered_map<std::string, std::string> &properties) {
        concord::Grid<uint8_t> grid(height, width, resolution_, true, shift_);
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

    void Raster::removeGrid(size_t index) {
        if (index < grid_layers_.size()) {
            grid_layers_.erase(grid_layers_.begin() + index);
        }
    }

    void Raster::addTerrainGrid(uint32_t width, uint32_t height, const std::string &name) {
        addGrid(width, height, name, "terrain");
    }

    void Raster::addOcclusionGrid(uint32_t width, uint32_t height, const std::string &name) {
        addGrid(width, height, name, "occlusion");
    }

    void Raster::addElevationGrid(uint32_t width, uint32_t height, const std::string &name) {
        addGrid(width, height, name, "elevation");
    }

    std::vector<GridLayer> Raster::getGridsByType(const std::string &type) const {
        std::vector<GridLayer> result;
        for (const auto &layer : grid_layers_) {
            if (layer.type == type) {
                result.push_back(layer);
            }
        }
        return result;
    }

    std::vector<GridLayer> Raster::filterByProperty(const std::string &key, const std::string &value) const {
        std::vector<GridLayer> result;
        for (const auto &layer : grid_layers_) {
            auto it = layer.properties.find(key);
            if (it != layer.properties.end() && it->second == value) {
                result.push_back(layer);
            }
        }
        return result;
    }

    std::vector<std::string> Raster::getGridNames() const {
        std::vector<std::string> names;
        for (const auto &layer : grid_layers_) {
            names.push_back(layer.name);
        }
        return names;
    }

    const concord::Datum &Raster::getDatum() const { return datum_; }
    void Raster::setDatum(const concord::Datum &datum) { datum_ = datum; }

    const concord::Pose &Raster::getShift() const { return shift_; }
    void Raster::setShift(const concord::Pose &shift) { shift_ = shift; }

    double Raster::getResolution() const { return resolution_; }
    void Raster::setResolution(double resolution) { resolution_ = resolution; }

    void Raster::setGlobalProperty(const std::string &key, const std::string &value) {
        for (auto &layer : grid_layers_) {
            layer.setGlobalProperty(key, value);
        }
    }

    std::string Raster::getGlobalProperty(const std::string &key, const std::string &default_value) const {
        if (!grid_layers_.empty()) {
            auto props = grid_layers_[0].getGlobalProperties();
            auto it = props.find(key);
            return (it != props.end()) ? it->second : default_value;
        }
        return default_value;
    }

    std::unordered_map<std::string, std::string> Raster::getGlobalProperties() const {
        return grid_layers_.empty() ? std::unordered_map<std::string, std::string>{}
                                    : grid_layers_[0].getGlobalProperties();
    }

    void Raster::removeGlobalProperty(const std::string &key) {
        std::hash<std::string> hasher;
        uint16_t tag = GLOBAL_PROPERTIES_BASE_TAG + (hasher(key) % 1000);
        for (auto &layer : grid_layers_) {
            layer.customTags.erase(tag);
        }
    }

} // namespace geotiv
