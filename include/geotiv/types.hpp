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
#include <pigment/pigment.hpp>

namespace dp = datapod;
namespace pg = pigment;

namespace geotiv {

    using std::uint16_t;

    /// Tag number ranges per TIFF 6.0 specification
    constexpr uint16_t TIFF_RESERVED_MAX = 32767;   // 0-32767 reserved for TIFF standard
    constexpr uint16_t PRIVATE_TAG_MIN = 32768;     // 32768-65535 for private/custom use
    constexpr uint16_t GEOTIV_RESERVED_MIN = 50000; // geotiv reserved range start
    constexpr uint16_t GEOTIV_RESERVED_MAX = 50999; // geotiv reserved range end
    constexpr uint16_t GLOBAL_PROPERTIES_BASE_TAG = 50100;

    /// Check if a tag number is safe for custom use
    /// Returns true if tag is in private range (>=32768) or geotiv reserved range (50000-50999)
    inline bool is_valid_custom_tag(uint16_t tag) {
        // Private range: 32768-65535
        if (tag >= PRIVATE_TAG_MIN) {
            return true;
        }
        // geotiv reserved range: 50000-50999 (subset of private range, but explicit)
        if (tag >= GEOTIV_RESERVED_MIN && tag <= GEOTIV_RESERVED_MAX) {
            return true;
        }
        return false;
    }

    /// Check if a tag number is in the TIFF reserved range
    inline bool is_reserved_tiff_tag(uint16_t tag) { return tag <= TIFF_RESERVED_MAX; }

    /// Validate a custom tag number, throwing if invalid
    inline void validate_custom_tag(uint16_t tag) {
        if (!is_valid_custom_tag(tag)) {
            throw std::runtime_error("Invalid custom tag number " + std::to_string(tag) +
                                     ": must be >= 32768 (private range). "
                                     "Tags 0-32767 are reserved for TIFF standard.");
        }
    }

    /// SampleFormat values as defined in TIFF spec (tag 339)
    enum class SampleFormat : uint16_t {
        UnsignedInt = 1, // Unsigned integer data (default)
        SignedInt = 2,   // Signed integer data (two's complement)
        Float = 3,       // IEEE floating point
        Undefined = 4    // Undefined/untyped data
    };

    /// PhotometricInterpretation values as defined in TIFF spec (tag 262)
    enum class PhotometricInterpretation : uint16_t {
        WhiteIsZero = 0, // Inverted grayscale (0 = white)
        BlackIsZero = 1, // Standard grayscale (0 = black)
        RGB = 2,         // RGB color image
        Palette = 3,     // Palette/indexed color (requires ColorMap)
        Mask = 4,        // Transparency mask
        CMYK = 5,        // CMYK color
        YCbCr = 6,       // YCbCr color (JPEG)
        CIELab = 8       // CIE L*a*b* color
    };

    /// RGBA color type using pigment library
    /// Stores R, G, B, A as uint8_t values (0-255)
    struct RGBA {
        uint8_t r = 0;
        uint8_t g = 0;
        uint8_t b = 0;
        uint8_t a = 255;

        RGBA() = default;
        RGBA(uint8_t r_, uint8_t g_, uint8_t b_, uint8_t a_ = 255) : r(r_), g(g_), b(b_), a(a_) {}

        /// Construct from pigment::RGB
        explicit RGBA(const pg::RGB &rgb, uint8_t alpha = 255) : r(rgb.r()), g(rgb.g()), b(rgb.b()), a(alpha) {}

        /// Convert to pigment::RGB
        pg::RGB to_rgb() const { return pg::RGB(r, g, b); }

        /// Check equality
        bool operator==(const RGBA &other) const {
            return r == other.r && g == other.g && b == other.b && a == other.a;
        }
        bool operator!=(const RGBA &other) const { return !(*this == other); }
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
                                     dp::Grid<double>,   // 64-bit float (high precision)
                                     dp::Grid<RGBA>      // RGBA color (RGB/RGBA images)
                                     >;

    /// Helper trait to extract element type from Grid<T>
    template <typename GridType> struct grid_element_type;
    template <typename T> struct grid_element_type<dp::Grid<T>> {
        using type = T;
    };
    template <typename GridType> using grid_element_type_t = typename grid_element_type<GridType>::type;

    /// Get bits per sample for a grid variant
    /// For RGBA grids, returns 8 (bits per channel)
    inline uint16_t get_bits_per_sample(const GridVariant &grid) {
        return std::visit(
            [](const auto &g) -> uint16_t {
                using GridType = std::decay_t<decltype(g)>;
                using T = grid_element_type_t<GridType>;
                if constexpr (std::is_same_v<T, RGBA>) {
                    return 8; // 8 bits per channel
                } else {
                    return static_cast<uint16_t>(sizeof(T) * 8);
                }
            },
            grid);
    }

    /// Get sample format for a grid variant
    inline SampleFormat get_sample_format(const GridVariant &grid) {
        return std::visit(
            [](const auto &g) -> SampleFormat {
                using GridType = std::decay_t<decltype(g)>;
                using T = grid_element_type_t<GridType>;
                if constexpr (std::is_same_v<T, RGBA>) {
                    return SampleFormat::UnsignedInt; // RGBA uses unsigned 8-bit per channel
                } else if constexpr (std::is_floating_point_v<T>) {
                    return SampleFormat::Float;
                } else if constexpr (std::is_signed_v<T>) {
                    return SampleFormat::SignedInt;
                } else {
                    return SampleFormat::UnsignedInt;
                }
            },
            grid);
    }

    /// Get samples per pixel for a grid variant
    /// Returns 4 for RGBA, 1 for scalar types
    inline uint16_t get_samples_per_pixel(const GridVariant &grid) {
        return std::visit(
            [](const auto &g) -> uint16_t {
                using GridType = std::decay_t<decltype(g)>;
                using T = grid_element_type_t<GridType>;
                if constexpr (std::is_same_v<T, RGBA>) {
                    return 4; // R, G, B, A
                } else {
                    return 1; // Single channel
                }
            },
            grid);
    }

    /// Get photometric interpretation for a grid variant
    inline PhotometricInterpretation get_photometric_interpretation(const GridVariant &grid) {
        return std::visit(
            [](const auto &g) -> PhotometricInterpretation {
                using GridType = std::decay_t<decltype(g)>;
                using T = grid_element_type_t<GridType>;
                if constexpr (std::is_same_v<T, RGBA>) {
                    return PhotometricInterpretation::RGB;
                } else {
                    return PhotometricInterpretation::BlackIsZero;
                }
            },
            grid);
    }

    /// Check if grid is a color (RGBA) grid
    inline bool is_color_grid(const GridVariant &grid) { return std::holds_alternative<dp::Grid<RGBA>>(grid); }

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

        /// Get photometric interpretation for this layer's grid
        inline PhotometricInterpretation photometricInterpretation() const {
            return get_photometric_interpretation(grid);
        }

        /// Check if this layer contains color (RGBA) data
        inline bool isColorLayer() const { return is_color_grid(grid); }

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

        /// Set a custom tag with validation
        /// Throws std::runtime_error if tag is in reserved TIFF range (0-32767)
        inline void setCustomTag(uint16_t tag, const std::vector<uint32_t> &values) {
            validate_custom_tag(tag);
            customTags[tag] = values;
        }

        /// Set a custom tag with a single value
        inline void setCustomTag(uint16_t tag, uint32_t value) { setCustomTag(tag, std::vector<uint32_t>{value}); }

        /// Get a custom tag value (returns empty vector if not found)
        inline std::vector<uint32_t> getCustomTag(uint16_t tag) const {
            auto it = customTags.find(tag);
            if (it != customTags.end()) {
                return it->second;
            }
            return {};
        }

        /// Check if a custom tag exists
        inline bool hasCustomTag(uint16_t tag) const { return customTags.find(tag) != customTags.end(); }

        inline void setGlobalProperty(const std::string &key, const std::string &value) {
            std::hash<std::string> hasher;
            uint16_t tag = GLOBAL_PROPERTIES_BASE_TAG + (hasher(key) % 1000);
            // Global properties are in geotiv reserved range, no validation needed
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
