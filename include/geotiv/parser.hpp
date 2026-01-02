#pragma once

#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "geotiv/types.hpp"
#include <concord/concord.hpp>
#include <datapod/datapod.hpp>

namespace dp = datapod;
namespace cc = concord;

namespace geotiv {

    namespace detail {
        struct TIFFEntry {
            uint16_t tag, type;
            uint32_t count, valueOffset;
        };

        inline uint16_t readLE16(std::ifstream &f) {
            uint16_t v = 0;
            f.read(reinterpret_cast<char *>(&v), sizeof(v));
            if (f.gcount() != sizeof(v)) {
                throw std::runtime_error("Failed to read LE16");
            }
            return v;
        }

        inline uint32_t readLE32(std::ifstream &f) {
            uint32_t v = 0;
            f.read(reinterpret_cast<char *>(&v), sizeof(v));
            if (f.gcount() != sizeof(v))
                throw std::runtime_error("Failed to read LE32");
            return v;
        }

        inline uint64_t readLE64(std::ifstream &f) {
            uint64_t v = 0;
            f.read(reinterpret_cast<char *>(&v), sizeof(v));
            if (f.gcount() != sizeof(v))
                throw std::runtime_error("Failed to read LE64");
            return v;
        }

        inline uint16_t readBE16(std::ifstream &f) {
            uint8_t bytes[2];
            f.read(reinterpret_cast<char *>(bytes), 2);
            if (f.gcount() != 2)
                throw std::runtime_error("Failed to read BE16");
            return (uint16_t(bytes[0]) << 8) | uint16_t(bytes[1]);
        }

        inline uint32_t readBE32(std::ifstream &f) {
            uint8_t bytes[4];
            f.read(reinterpret_cast<char *>(bytes), 4);
            if (f.gcount() != 4)
                throw std::runtime_error("Failed to read BE32");
            return (uint32_t(bytes[0]) << 24) | (uint32_t(bytes[1]) << 16) | (uint32_t(bytes[2]) << 8) |
                   uint32_t(bytes[3]);
        }

        inline uint64_t readBE64(std::ifstream &f) {
            uint8_t bytes[8];
            f.read(reinterpret_cast<char *>(bytes), 8);
            if (f.gcount() != 8)
                throw std::runtime_error("Failed to read BE64");
            uint64_t result = 0;
            for (int i = 0; i < 8; ++i) {
                result = (result << 8) | uint64_t(bytes[i]);
            }
            return result;
        }

        inline std::string readString(std::ifstream &f, uint32_t offset, uint32_t count) {
            if (count == 0)
                return "";
            std::vector<char> buf(count);
            f.seekg(offset, std::ios::beg);
            f.read(buf.data(), count);
            if (f.gcount() != static_cast<std::streamsize>(count))
                throw std::runtime_error("Failed to read string data");

            // Ensure null termination
            if (buf.back() != '\0') {
                buf.push_back('\0');
            }
            return std::string(buf.data());
        }
    } // namespace detail

    namespace fs = std::filesystem;

    inline RasterCollection ReadRasterCollection(const fs::path &file) {
        std::ifstream f(file, std::ios::binary);
        if (!f)
            throw std::runtime_error("Cannot open \"" + file.string() + "\"");

        // 1) Header
        char bom[2];
        f.read(bom, 2);
        bool little = (bom[0] == 'I' && bom[1] == 'I');
        if (!little && !(bom[0] == 'M' && bom[1] == 'M'))
            throw std::runtime_error("Bad TIFF byte-order");

        auto read16 = little ? detail::readLE16 : detail::readBE16;
        auto read32 = little ? detail::readLE32 : detail::readBE32;
        auto read64 = little ? detail::readLE64 : detail::readBE64;

        if (read16(f) != 42)
            throw std::runtime_error("Bad TIFF magic");
        uint32_t nextIFD = read32(f);

        RasterCollection rc;

        // 2) Loop IFDs
        bool firstIFD = true;
        while (nextIFD) {
            f.seekg(nextIFD, std::ios::beg);
            uint16_t nEnt = read16(f);
            std::map<uint16_t, detail::TIFFEntry> E;
            for (int i = 0; i < nEnt; ++i) {
                detail::TIFFEntry e;
                e.tag = read16(f);
                e.type = read16(f);
                e.count = read32(f);
                e.valueOffset = read32(f);
                E[e.tag] = e;
            }
            uint32_t currentIFDOffset = nextIFD;
            nextIFD = read32(f);

            // TIFF type constants
            constexpr uint16_t TIFF_BYTE = 1;       // 8-bit unsigned
            constexpr uint16_t TIFF_ASCII = 2;      // 8-bit ASCII
            constexpr uint16_t TIFF_SHORT = 3;      // 16-bit unsigned
            constexpr uint16_t TIFF_LONG = 4;       // 32-bit unsigned
            constexpr uint16_t TIFF_RATIONAL = 5;   // Two LONGs (num/denom)
            constexpr uint16_t TIFF_SBYTE = 6;      // 8-bit signed
            constexpr uint16_t TIFF_UNDEFINED = 7;  // 8-bit untyped
            constexpr uint16_t TIFF_SSHORT = 8;     // 16-bit signed
            constexpr uint16_t TIFF_SLONG = 9;      // 32-bit signed
            constexpr uint16_t TIFF_SRATIONAL = 10; // Two SLONGs
            constexpr uint16_t TIFF_FLOAT = 11;     // 32-bit IEEE float
            constexpr uint16_t TIFF_DOUBLE = 12;    // 64-bit IEEE double

            // Helper: get single unsigned integer value from tag
            // Supports BYTE, SHORT, LONG types
            auto getUInt = [&](uint16_t tag) -> uint32_t {
                auto it = E.find(tag);
                if (it == E.end())
                    return 0;
                auto &e = it->second;

                if (e.type == TIFF_BYTE || e.type == TIFF_UNDEFINED) {
                    // BYTE: up to 4 values fit in valueOffset field
                    if (e.count >= 1) {
                        return little ? (e.valueOffset & 0xFF) : ((e.valueOffset >> 24) & 0xFF);
                    }
                } else if (e.type == TIFF_SHORT) {
                    if (e.count == 1) {
                        return little ? (e.valueOffset & 0xFFFF) : ((e.valueOffset >> 16) & 0xFFFF);
                    } else {
                        // Multiple shorts - read from offset
                        f.seekg(e.valueOffset, std::ios::beg);
                        return read16(f);
                    }
                } else if (e.type == TIFF_LONG) {
                    return e.valueOffset;
                } else if (e.type == TIFF_SBYTE) {
                    // Signed byte - return as unsigned for compatibility
                    if (e.count >= 1) {
                        return little ? (e.valueOffset & 0xFF) : ((e.valueOffset >> 24) & 0xFF);
                    }
                } else if (e.type == TIFF_SSHORT) {
                    if (e.count == 1) {
                        return little ? (e.valueOffset & 0xFFFF) : ((e.valueOffset >> 16) & 0xFFFF);
                    }
                } else if (e.type == TIFF_SLONG) {
                    return e.valueOffset;
                }
                return 0;
            };

            // Helper: read multiple unsigned integer values from tag
            // Supports BYTE, SHORT, LONG types
            auto readUInts = [&](uint16_t tag) -> std::vector<uint32_t> {
                auto it = E.find(tag);
                if (it == E.end())
                    return {};
                auto &e = it->second;

                std::vector<uint32_t> out;

                if (e.type == TIFF_BYTE || e.type == TIFF_UNDEFINED || e.type == TIFF_SBYTE) {
                    // BYTE: up to 4 values fit in valueOffset field
                    if (e.count <= 4) {
                        for (uint32_t i = 0; i < e.count; ++i) {
                            uint8_t val;
                            if (little) {
                                val = (e.valueOffset >> (8 * i)) & 0xFF;
                            } else {
                                val = (e.valueOffset >> (8 * (3 - i))) & 0xFF;
                            }
                            out.push_back(val);
                        }
                    } else {
                        f.seekg(e.valueOffset, std::ios::beg);
                        for (uint32_t i = 0; i < e.count; ++i) {
                            uint8_t val;
                            f.read(reinterpret_cast<char *>(&val), 1);
                            out.push_back(val);
                        }
                    }
                } else if (e.type == TIFF_SHORT || e.type == TIFF_SSHORT) {
                    if (e.count == 1) {
                        uint16_t val = little ? (e.valueOffset & 0xFFFF) : ((e.valueOffset >> 16) & 0xFFFF);
                        out.push_back(val);
                    } else if (e.count == 2) {
                        uint16_t val1 = little ? (e.valueOffset & 0xFFFF) : ((e.valueOffset >> 16) & 0xFFFF);
                        uint16_t val2 = little ? ((e.valueOffset >> 16) & 0xFFFF) : (e.valueOffset & 0xFFFF);
                        out.push_back(val1);
                        out.push_back(val2);
                    } else {
                        f.seekg(e.valueOffset, std::ios::beg);
                        for (uint32_t i = 0; i < e.count; ++i)
                            out.push_back(read16(f));
                    }
                } else if (e.type == TIFF_LONG || e.type == TIFF_SLONG) {
                    if (e.count == 1) {
                        out.push_back(e.valueOffset);
                    } else {
                        f.seekg(e.valueOffset, std::ios::beg);
                        for (uint32_t i = 0; i < e.count; ++i)
                            out.push_back(read32(f));
                    }
                }
                return out;
            };

            // Helper: read RATIONAL values (pairs of LONGs as numerator/denominator)
            // Returns as doubles for convenience
            auto readRationals = [&](uint16_t tag) -> std::vector<double> {
                auto it = E.find(tag);
                if (it == E.end())
                    return {};
                auto &e = it->second;

                std::vector<double> out;
                if (e.type != TIFF_RATIONAL && e.type != TIFF_SRATIONAL)
                    return out;

                f.seekg(e.valueOffset, std::ios::beg);
                for (uint32_t i = 0; i < e.count; ++i) {
                    uint32_t num = read32(f);
                    uint32_t denom = read32(f);
                    if (denom == 0) {
                        out.push_back(0.0); // Avoid division by zero
                    } else if (e.type == TIFF_SRATIONAL) {
                        out.push_back(static_cast<double>(static_cast<int32_t>(num)) /
                                      static_cast<double>(static_cast<int32_t>(denom)));
                    } else {
                        out.push_back(static_cast<double>(num) / static_cast<double>(denom));
                    }
                }
                return out;
            };

            // Helper: read FLOAT values (32-bit IEEE)
            auto readFloats = [&](uint16_t tag) -> std::vector<float> {
                auto it = E.find(tag);
                if (it == E.end())
                    return {};
                auto &e = it->second;

                std::vector<float> out;
                if (e.type != TIFF_FLOAT)
                    return out;

                if (e.count == 1) {
                    // Single float fits in valueOffset
                    float fval;
                    std::memcpy(&fval, &e.valueOffset, sizeof(fval));
                    out.push_back(fval);
                } else {
                    f.seekg(e.valueOffset, std::ios::beg);
                    for (uint32_t i = 0; i < e.count; ++i) {
                        uint32_t bits = read32(f);
                        float fval;
                        std::memcpy(&fval, &bits, sizeof(fval));
                        out.push_back(fval);
                    }
                }
                return out;
            };

            // Helper: read DOUBLE values (64-bit IEEE)
            auto readDoubles = [&](uint16_t tag) -> std::vector<double> {
                auto it = E.find(tag);
                if (it == E.end())
                    return {};
                auto &e = it->second;

                std::vector<double> out;

                // Also accept FLOAT type and convert to double
                if (e.type == TIFF_FLOAT) {
                    auto floats = readFloats(tag);
                    for (float fv : floats) {
                        out.push_back(static_cast<double>(fv));
                    }
                    return out;
                }

                // Also accept RATIONAL type
                if (e.type == TIFF_RATIONAL || e.type == TIFF_SRATIONAL) {
                    return readRationals(tag);
                }

                if (e.type != TIFF_DOUBLE)
                    return out;

                f.seekg(e.valueOffset, std::ios::beg);
                for (uint32_t i = 0; i < e.count; ++i) {
                    uint64_t bits = read64(f);
                    double d;
                    std::memcpy(&d, &bits, sizeof(d));
                    out.push_back(d);
                }
                return out;
            };

            // build Layer - validate required tags
            Layer L;
            L.ifdOffset = currentIFDOffset;
            L.width = getUInt(256);  // ImageWidth
            L.height = getUInt(257); // ImageLength

            if (L.width == 0 || L.height == 0)
                throw std::runtime_error("Invalid or missing image dimensions");

            L.samplesPerPixel = getUInt(277);
            if (L.samplesPerPixel == 0)
                L.samplesPerPixel = 1; // default

            L.planarConfig = getUInt(284);
            if (L.planarConfig == 0)
                L.planarConfig = 1; // default: chunky

            uint32_t bitsPerSample = getUInt(258);
            if (bitsPerSample == 0)
                bitsPerSample = 8; // default to 8-bit

            // Read SampleFormat tag (339) - defaults to 1 (unsigned int) if not present
            uint32_t sampleFormatValue = getUInt(339);
            if (sampleFormatValue == 0)
                sampleFormatValue = 1; // default: unsigned integer
            SampleFormat sampleFormat = static_cast<SampleFormat>(sampleFormatValue);

            // Read PhotometricInterpretation tag (262) - defaults to 1 (BlackIsZero) if not present
            uint32_t photometricValue = getUInt(262);
            PhotometricInterpretation photometric = static_cast<PhotometricInterpretation>(photometricValue);

            // Validate supported bit depths
            if (bitsPerSample != 8 && bitsPerSample != 16 && bitsPerSample != 32 && bitsPerSample != 64) {
                throw std::runtime_error("Unsupported bits per sample: " + std::to_string(bitsPerSample) +
                                         ". Supported: 8, 16, 32, 64");
            }

            // Compression and Predictor tags are ignored (not supported)

            L.stripOffsets = readUInts(273);    // StripOffsets
            L.stripByteCounts = readUInts(279); // StripByteCounts

            if (L.stripOffsets.empty() || L.stripByteCounts.empty())
                throw std::runtime_error("Missing strip data");
            if (L.stripOffsets.size() != L.stripByteCounts.size())
                throw std::runtime_error("Mismatched strip arrays");

            // Read all strips and combine pixel data
            size_t totalBytes = 0;
            for (auto count : L.stripByteCounts) {
                totalBytes += count;
            }

            size_t bytesPerSample = bitsPerSample / 8;
            size_t expectedBytes = size_t(L.width) * L.height * L.samplesPerPixel * bytesPerSample;

            // Verify size matches (only uncompressed data supported)
            if (totalBytes != expectedBytes) {
                throw std::runtime_error("Strip byte count mismatch: expected " + std::to_string(expectedBytes) +
                                         ", got " + std::to_string(totalBytes));
            }

            // Read strip data
            std::vector<uint8_t> pix(totalBytes);
            size_t offset = 0;

            for (size_t i = 0; i < L.stripOffsets.size(); ++i) {
                f.seekg(L.stripOffsets[i], std::ios::beg);
                f.read(reinterpret_cast<char *>(pix.data() + offset), L.stripByteCounts[i]);
                if (f.gcount() != static_cast<std::streamsize>(L.stripByteCounts[i]))
                    throw std::runtime_error("Failed to read strip data");
                offset += L.stripByteCounts[i];
            }

            // Parse geotags for each IFD independently (always WGS84)
            dp::Geo layerDatum;           // Will be set from ImageDescription or use a valid default
            dp::Pose layerShift{};        // default (identity quaternion)
            double layerResolution = 1.0; // default
            std::string layerDescription;
            bool datumFromDescription = false;

            // Parse ImageDescription for CRS/DATUM/HEADING
            auto itD = E.find(270);
            if (itD != E.end() && itD->second.type == 2) {
                layerDescription = detail::readString(f, itD->second.valueOffset, itD->second.count);
                std::istringstream ss(layerDescription);
                std::string tok;

                while (ss >> tok) {
                    if (tok == "CRS") {
                        std::string s;
                        if (ss >> s) {
                            // All CRS are WGS84 - ignore the parsed value
                        }
                    } else if (tok == "DATUM") {
                        if (ss >> layerDatum.latitude >> layerDatum.longitude >> layerDatum.altitude) {
                            datumFromDescription = true;
                        }
                    } else if (tok == "SHIFT") {
                        double x, y, z, yaw;
                        if (ss >> x >> y >> z >> yaw) {
                            layerShift = dp::Pose{dp::Point{x, y, z}, dp::Quaternion::from_euler(0, 0, yaw)};
                        }
                    }
                }
            }

            // If we didn't get a valid datum from ImageDescription, use a default that passes is_set()
            if (!datumFromDescription) {
                layerDatum = dp::Geo{0.001, 0.001, 1.0}; // Valid minimal coordinates
            }

            // Try ModelTransformationTag (34264) first for rotated grids
            // Then fall back to ModelPixelScaleTag (33550) for non-rotated grids
            auto transform = readDoubles(34264);
            if (transform.size() >= 16) {
                // ModelTransformationTag: 4x4 affine matrix (row-major)
                // [a b 0 d; e f 0 h; 0 0 1 0; 0 0 0 1]
                // a = scale_x * cos(θ), b = -scale_y * sin(θ)
                // e = scale_x * sin(θ), f = scale_y * cos(θ)
                double a = transform[0];
                double b = transform[1];
                double e = transform[4];
                double f = transform[5];

                // Extract scale and rotation from the matrix
                // scale_x = sqrt(a² + e²) gives the X scale in degrees
                double scale_x_deg = std::sqrt(a * a + e * e);

                // Extract rotation angle: θ = atan2(e, a)
                double yaw = std::atan2(e, a);
                layerShift.rotation = dp::Quaternion::from_euler(0, 0, yaw);

                // Convert X scale from degrees to meters using precise concord conversion
                // Use the X scale (longitude) since that's what we use for non-rotated grids too
                cc::earth::WGS center_wgs{layerDatum.latitude, layerDatum.longitude, layerDatum.altitude};
                cc::earth::WGS east_point_wgs{layerDatum.latitude, layerDatum.longitude + scale_x_deg,
                                              layerDatum.altitude};

                cc::frame::ENU center_enu = cc::frame::to_enu(layerDatum, center_wgs);
                cc::frame::ENU east_point_enu = cc::frame::to_enu(layerDatum, east_point_wgs);

                layerResolution = east_point_enu.x() - center_enu.x();
            } else {
                // ModelPixelScale → resolution for this IFD (non-rotated grids)
                auto scales = readDoubles(33550);
                if (scales.size() >= 2) {
                    // PixelScale is stored in degrees, use precise concord conversions to get back meters
                    double resolution_deg_lon = scales[0]; // X scale in degrees
                    // scales[1] is Y scale (negative for standard GeoTIFF), but we use X scale for resolution

                    // Use precise concord conversion: create two WGS points and convert to ENU to get distance
                    cc::earth::WGS center_wgs{layerDatum.latitude, layerDatum.longitude, layerDatum.altitude};
                    cc::earth::WGS east_point_wgs{layerDatum.latitude, layerDatum.longitude + resolution_deg_lon,
                                                  layerDatum.altitude};

                    cc::frame::ENU center_enu = cc::frame::to_enu(layerDatum, center_wgs);
                    cc::frame::ENU east_point_enu = cc::frame::to_enu(layerDatum, east_point_wgs);

                    // Calculate precise meter distance (use absolute value since Y scale is negative)
                    layerResolution = east_point_enu.x() - center_enu.x();
                }
            }

            if (layerResolution <= 0) {
                throw std::runtime_error("Invalid pixel scale: " + std::to_string(layerResolution));
            }

            // Set layer-specific metadata (always WGS84)
            L.datum = layerDatum;
            L.shift = layerShift;
            L.resolution = layerResolution;
            L.imageDescription = layerDescription;

            // Read custom tags (tag numbers 50000 and above are typically custom)
            for (const auto &[tag, entry] : E) {
                if (tag >= 50000) {
                    L.customTags[tag] = readUInts(tag);
                }
            }

            // Set collection defaults from first IFD if not set
            if (firstIFD) {
                firstIFD = false;
                rc.datum = layerDatum;
                rc.shift = layerShift;
                rc.resolution = layerResolution;
            }

            // Build geo-grid using layer-specific resolution and datum
            if (!L.datum.is_set()) {
                throw std::runtime_error("Datum not properly initialized for layer");
            }

            // Use the shift directly - it's already in ENU space
            dp::Pose shift = L.shift;

            // Helper lambda to read a value from the pixel buffer in little-endian format
            auto readLE = [&pix, little](size_t offset, size_t bytes) -> uint64_t {
                uint64_t value = 0;
                if (little) {
                    for (size_t i = 0; i < bytes; ++i) {
                        value |= static_cast<uint64_t>(pix[offset + i]) << (8 * i);
                    }
                } else {
                    for (size_t i = 0; i < bytes; ++i) {
                        value |= static_cast<uint64_t>(pix[offset + i]) << (8 * (bytes - 1 - i));
                    }
                }
                return value;
            };

            // Create and fill grid based on data type
            // Lambda to fill grid from pixel data
            auto fillGrid = [&](auto &grid, size_t bytesPerSample) {
                using GridType = std::decay_t<decltype(grid)>;
                using T = grid_element_type_t<GridType>;
                if (L.planarConfig == 1) { // Chunky format
                    size_t idx = 0;
                    for (int32_t r = L.height - 1; r >= 0; --r) {
                        for (uint32_t c = 0; c < L.width; ++c) {
                            uint64_t rawValue = readLE(idx, bytesPerSample);
                            if constexpr (std::is_floating_point_v<T>) {
                                if constexpr (sizeof(T) == 4) {
                                    uint32_t bits = static_cast<uint32_t>(rawValue);
                                    float fval;
                                    std::memcpy(&fval, &bits, sizeof(fval));
                                    grid(r, c) = fval;
                                } else {
                                    double dval;
                                    std::memcpy(&dval, &rawValue, sizeof(dval));
                                    grid(r, c) = dval;
                                }
                            } else {
                                grid(r, c) = static_cast<T>(rawValue);
                            }
                            idx += bytesPerSample * L.samplesPerPixel;
                        }
                    }
                } else { // Planar format
                    size_t idx = 0;
                    for (int32_t r = L.height - 1; r >= 0; --r) {
                        for (uint32_t c = 0; c < L.width; ++c) {
                            uint64_t rawValue = readLE(idx, bytesPerSample);
                            if constexpr (std::is_floating_point_v<T>) {
                                if constexpr (sizeof(T) == 4) {
                                    uint32_t bits = static_cast<uint32_t>(rawValue);
                                    float fval;
                                    std::memcpy(&fval, &bits, sizeof(fval));
                                    grid(r, c) = fval;
                                } else {
                                    double dval;
                                    std::memcpy(&dval, &rawValue, sizeof(dval));
                                    grid(r, c) = dval;
                                }
                            } else {
                                grid(r, c) = static_cast<T>(rawValue);
                            }
                            idx += bytesPerSample;
                        }
                    }
                }
            };

            // Create appropriate grid type based on bitsPerSample, sampleFormat, and photometric
            // Check for RGB/RGBA images first (PhotometricInterpretation = 2)
            if (photometric == PhotometricInterpretation::RGB && bitsPerSample == 8) {
                // RGB or RGBA image - create RGBA grid
                auto grid = dp::make_grid<RGBA>(L.height, L.width, L.resolution, true, shift, RGBA{});

                if (L.planarConfig == 1) { // Chunky format (RGBRGB... or RGBARGBA...)
                    size_t idx = 0;
                    size_t bytesPerPixel = L.samplesPerPixel; // 3 for RGB, 4 for RGBA
                    for (int32_t r = L.height - 1; r >= 0; --r) {
                        for (uint32_t c = 0; c < L.width; ++c) {
                            RGBA pixel;
                            pixel.r = pix[idx];
                            pixel.g = pix[idx + 1];
                            pixel.b = pix[idx + 2];
                            pixel.a = (L.samplesPerPixel >= 4) ? pix[idx + 3] : 255;
                            grid(r, c) = pixel;
                            idx += bytesPerPixel;
                        }
                    }
                } else { // Planar format (RRR...GGG...BBB...AAA...)
                    size_t planeSize = size_t(L.width) * L.height;
                    for (int32_t r = L.height - 1; r >= 0; --r) {
                        for (uint32_t c = 0; c < L.width; ++c) {
                            size_t pixelIdx = size_t(L.height - 1 - r) * L.width + c;
                            RGBA pixel;
                            pixel.r = pix[pixelIdx];
                            pixel.g = pix[planeSize + pixelIdx];
                            pixel.b = pix[2 * planeSize + pixelIdx];
                            pixel.a = (L.samplesPerPixel >= 4) ? pix[3 * planeSize + pixelIdx] : 255;
                            grid(r, c) = pixel;
                        }
                    }
                }
                L.grid = std::move(grid);
            } else if (bitsPerSample == 8) {
                if (sampleFormat == SampleFormat::SignedInt) {
                    auto grid = dp::make_grid<int8_t>(L.height, L.width, L.resolution, true, shift, int8_t{0});
                    fillGrid(grid, 1);
                    L.grid = std::move(grid);
                } else {
                    auto grid = dp::make_grid<uint8_t>(L.height, L.width, L.resolution, true, shift, uint8_t{0});
                    fillGrid(grid, 1);
                    L.grid = std::move(grid);
                }
            } else if (bitsPerSample == 16) {
                if (sampleFormat == SampleFormat::SignedInt) {
                    auto grid = dp::make_grid<int16_t>(L.height, L.width, L.resolution, true, shift, int16_t{0});
                    fillGrid(grid, 2);
                    L.grid = std::move(grid);
                } else {
                    auto grid = dp::make_grid<uint16_t>(L.height, L.width, L.resolution, true, shift, uint16_t{0});
                    fillGrid(grid, 2);
                    L.grid = std::move(grid);
                }
            } else if (bitsPerSample == 32) {
                if (sampleFormat == SampleFormat::Float) {
                    auto grid = dp::make_grid<float>(L.height, L.width, L.resolution, true, shift, 0.0f);
                    fillGrid(grid, 4);
                    L.grid = std::move(grid);
                } else if (sampleFormat == SampleFormat::SignedInt) {
                    auto grid = dp::make_grid<int32_t>(L.height, L.width, L.resolution, true, shift, int32_t{0});
                    fillGrid(grid, 4);
                    L.grid = std::move(grid);
                } else {
                    auto grid = dp::make_grid<uint32_t>(L.height, L.width, L.resolution, true, shift, uint32_t{0});
                    fillGrid(grid, 4);
                    L.grid = std::move(grid);
                }
            } else if (bitsPerSample == 64) {
                if (sampleFormat == SampleFormat::Float) {
                    auto grid = dp::make_grid<double>(L.height, L.width, L.resolution, true, shift, 0.0);
                    fillGrid(grid, 8);
                    L.grid = std::move(grid);
                } else {
                    throw std::runtime_error("64-bit integer grids not supported");
                }
            }
            rc.layers.emplace_back(std::move(L));
        }

        if (rc.layers.empty()) {
            throw std::runtime_error("No valid IFDs found in TIFF file");
        }

        return rc;
    }

    inline std::ostream &operator<<(std::ostream &os, RasterCollection const &rc) {
        os << "GeoTIFF RasterCollection\n"
           << " CRS:        WGS84\n"
           << " DATUM:      " << rc.datum.latitude << ", " << rc.datum.longitude << ", " << rc.datum.altitude << "\n"
           << " SHIFT:      " << rc.shift.point.x << ", " << rc.shift.point.y << ", " << rc.shift.point.z
           << " (yaw=" << rc.shift.rotation.to_euler().yaw << ")\n"
           << " RESOLUTION: " << rc.resolution << " (map units per pixel)\n"
           << " Layers:     " << rc.layers.size() << "\n";
        for (auto const &L : rc.layers) {
            os << "  IFD@0x" << std::hex << L.ifdOffset << std::dec << " → " << L.width << "×" << L.height
               << ", SPP=" << L.samplesPerPixel << ", PC=" << L.planarConfig << "\n";
        }
        return os;
    }

} // namespace geotiv
