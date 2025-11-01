#pragma once

#include <filesystem>
#include <functional>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "concord/concord.hpp"

namespace geotiv {

    using std::uint16_t;

    constexpr uint16_t GLOBAL_PROPERTIES_BASE_TAG = 50100;

    std::vector<uint32_t> stringToAsciiTag(const std::string &str);
    std::string asciiTagToString(const std::vector<uint32_t> &data);

    struct Layer {
        uint32_t ifdOffset = 0;

        uint32_t width = 0;
        uint32_t height = 0;
        uint32_t samplesPerPixel = 0;
        uint32_t planarConfig = 0;

        std::vector<uint32_t> stripOffsets;
        std::vector<uint32_t> stripByteCounts;

        concord::Datum datum;
        concord::Pose shift;
        double resolution = 1.0;

        std::string imageDescription;
        std::map<uint16_t, std::vector<uint32_t>> customTags;

        concord::Grid<uint8_t> grid;

        void setGlobalProperty(const std::string &key, const std::string &value);
        std::unordered_map<std::string, std::string> getGlobalProperties() const;
    };

    struct RasterCollection {
        std::vector<Layer> layers;

        concord::Datum datum;
        concord::Pose shift;
        double resolution;

        std::unordered_map<std::string, std::string> getGlobalPropertiesFromFirstLayer() const;
        void setGlobalPropertiesOnAllLayers(const std::unordered_map<std::string, std::string> &props);
    };

} // namespace geotiv
