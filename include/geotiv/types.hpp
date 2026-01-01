#pragma once

#include <cstdint>
#include <filesystem>
#include <functional>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

#include <datapod/datapod.hpp>

namespace dp = datapod;

namespace geotiv {

    using std::uint16_t;

    constexpr uint16_t GLOBAL_PROPERTIES_BASE_TAG = 50100;

    /// SampleFormat values as defined in TIFF spec (tag 339)
    enum class SampleFormat : uint16_t {
        UnsignedInt = 1, // Unsigned integer data (default)
        SignedInt = 2,   // Signed integer data (two's complement)
        Float = 3,       // IEEE floating point
        Undefined = 4    // Undefined/untyped data
    };

    /// Variant type supporting all grid data types
    /// This allows a Layer to hold grids of different numeric types
    using GridVariant = std::variant<dp::Grid<uint8_t>,  // 8-bit unsigned (default)
                                     dp::Grid<int8_t>,   // 8-bit signed
                                     dp::Grid<uint16_t>, // 16-bit unsigned (DEMs, elevation)
                                     dp::Grid<int16_t>,  // 16-bit signed (temperature, relative)
                                     dp::Grid<uint32_t>, // 32-bit unsigned (large counts)
                                     dp::Grid<int32_t>,  // 32-bit signed
                                     dp::Grid<float>,    // 32-bit float (scientific data)
                                     dp::Grid<double>    // 64-bit float (high precision)
                                     >;

    /// Helper trait to extract element type from Grid<T>
    template <typename GridType> struct grid_element_type;
    template <typename T> struct grid_element_type<dp::Grid<T>> {
        using type = T;
    };
    template <typename GridType> using grid_element_type_t = typename grid_element_type<GridType>::type;

    /// Get bits per sample for a grid variant
    inline uint16_t get_bits_per_sample(const GridVariant &grid) {
        return std::visit(
            [](const auto &g) -> uint16_t {
                using GridType = std::decay_t<decltype(g)>;
                using T = grid_element_type_t<GridType>;
                return static_cast<uint16_t>(sizeof(T) * 8);
            },
            grid);
    }

    /// Get sample format for a grid variant
    inline SampleFormat get_sample_format(const GridVariant &grid) {
        return std::visit(
            [](const auto &g) -> SampleFormat {
                using GridType = std::decay_t<decltype(g)>;
                using T = grid_element_type_t<GridType>;
                if constexpr (std::is_floating_point_v<T>) {
                    return SampleFormat::Float;
                } else if constexpr (std::is_signed_v<T>) {
                    return SampleFormat::SignedInt;
                } else {
                    return SampleFormat::UnsignedInt;
                }
            },
            grid);
    }

    /// Get grid dimensions (rows, cols)
    inline std::pair<size_t, size_t> get_grid_dimensions(const GridVariant &grid) {
        return std::visit([](const auto &g) -> std::pair<size_t, size_t> { return {g.rows, g.cols}; }, grid);
    }

    /// Get grid resolution
    inline double get_grid_resolution(const GridVariant &grid) {
        return std::visit([](const auto &g) -> double { return g.resolution; }, grid);
    }

    /// Get grid pose
    inline dp::Pose get_grid_pose(const GridVariant &grid) {
        return std::visit([](const auto &g) -> dp::Pose { return g.pose; }, grid);
    }

    /// Check if grid holds a specific type
    template <typename T> inline bool holds_grid_type(const GridVariant &grid) {
        return std::holds_alternative<dp::Grid<T>>(grid);
    }

    /// Get grid as specific type (throws if wrong type)
    template <typename T> inline dp::Grid<T> &get_grid(GridVariant &grid) { return std::get<dp::Grid<T>>(grid); }

    /// Get grid as specific type (const version)
    template <typename T> inline const dp::Grid<T> &get_grid(const GridVariant &grid) {
        return std::get<dp::Grid<T>>(grid);
    }

    /// Try to get grid as specific type (returns nullptr if wrong type)
    template <typename T> inline dp::Grid<T> *get_grid_if(GridVariant &grid) { return std::get_if<dp::Grid<T>>(&grid); }

    /// Try to get grid as specific type (const version)
    template <typename T> inline const dp::Grid<T> *get_grid_if(const GridVariant &grid) {
        return std::get_if<dp::Grid<T>>(&grid);
    }

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

        /// Grid data - supports multiple numeric types via variant
        GridVariant grid;

        /// Get bits per sample for this layer's grid
        inline uint16_t bitsPerSample() const { return get_bits_per_sample(grid); }

        /// Get sample format for this layer's grid
        inline SampleFormat sampleFormat() const { return get_sample_format(grid); }

        /// Check if grid holds a specific type
        template <typename T> inline bool holdsType() const { return holds_grid_type<T>(grid); }

        /// Get grid as specific type (throws std::bad_variant_access if wrong type)
        template <typename T> inline dp::Grid<T> &gridAs() { return get_grid<T>(grid); }

        /// Get grid as specific type (const version)
        template <typename T> inline const dp::Grid<T> &gridAs() const { return get_grid<T>(grid); }

        /// Try to get grid as specific type (returns nullptr if wrong type)
        template <typename T> inline dp::Grid<T> *gridIf() { return get_grid_if<T>(grid); }

        /// Try to get grid as specific type (const version)
        template <typename T> inline const dp::Grid<T> *gridIf() const { return get_grid_if<T>(grid); }

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
