#pragma once

#include <filesystem>
#include <functional>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include <datapod/datapod.hpp>

namespace dp = datapod;

namespace geotiv {

    using std::uint16_t;

    constexpr uint16_t GLOBAL_PROPERTIES_BASE_TAG = 50100;

    inline std::vector<uint32_t> stringToAsciiTag(const std::string &str) {
        std::string padded = str + '\0';
        while (padded.size() % 4 != 0) {
            padded += '\0';
        }

        std::vector<uint32_t> result;
        result.reserve(padded.size() / 4);

        for (size_t i = 0; i < padded.size(); i += 4) {
            uint32_t value = 0;
            for (int j = 0; j < 4; ++j) {
                value |= (static_cast<uint32_t>(padded[i + j]) << (j * 8));
            }
            result.push_back(value);
        }
        return result;
    }

    inline std::string asciiTagToString(const std::vector<uint32_t> &data) {
        std::string result;
        result.reserve(data.size() * 4);

        for (uint32_t value : data) {
            for (int j = 0; j < 4; ++j) {
                char c = static_cast<char>((value >> (j * 8)) & 0xFF);
                if (c == '\0') {
                    return result;
                }
                result += c;
            }
        }
        return result;
    }

    struct Layer {
        uint32_t ifdOffset = 0;

        uint32_t width = 0;
        uint32_t height = 0;
        uint32_t samplesPerPixel = 0;
        uint32_t planarConfig = 0;

        std::vector<uint32_t> stripOffsets;
        std::vector<uint32_t> stripByteCounts;

        dp::Geo datum;
        dp::Pose shift;
        double resolution = 1.0;

        std::string imageDescription;
        std::map<uint16_t, std::vector<uint32_t>> customTags;

        dp::Grid<uint8_t> grid;

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

    struct RasterCollection {
        std::vector<Layer> layers;

        dp::Geo datum;
        dp::Pose shift;
        double resolution;

        inline std::unordered_map<std::string, std::string> getGlobalPropertiesFromFirstLayer() const {
            if (!layers.empty()) {
                return layers[0].getGlobalProperties();
            }
            return {};
        }

        inline void setGlobalPropertiesOnAllLayers(const std::unordered_map<std::string, std::string> &props) {
            for (auto &layer : layers) {
                for (const auto &[key, value] : props) {
                    layer.setGlobalProperty(key, value);
                }
            }
        }
    };

} // namespace geotiv
