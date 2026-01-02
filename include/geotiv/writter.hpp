#pragma once

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "geotiv/types.hpp"
#include <concord/concord.hpp>
#include <datapod/datapod.hpp>

namespace dp = datapod;
namespace cc = concord;

namespace geotiv {

    namespace fs = std::filesystem;

    /// TIFF IFD entry structure for building sorted tag lists
    struct IfdEntry {
        uint16_t tag;
        uint16_t type;
        uint32_t count;
        uint32_t value_or_offset;

        bool operator<(const IfdEntry &other) const { return tag < other.tag; }
    };

    /// Check if a layer has significant rotation (non-zero yaw)
    /// Returns true if rotation should be encoded using ModelTransformationTag
    inline bool has_rotation(const Layer &layer) {
        auto euler = layer.shift.rotation.to_euler();
        // Consider rotation significant if yaw is more than ~0.01 degrees
        constexpr double ROTATION_THRESHOLD = 0.0002; // ~0.01 degrees in radians
        return std::abs(euler.yaw) > ROTATION_THRESHOLD;
    }

    /// Helper to write a value in little-endian format to a buffer
    template <typename T> inline void write_le(std::vector<uint8_t> &buf, size_t &pos, T value) {
        for (size_t i = 0; i < sizeof(T); ++i) {
            buf[pos++] = static_cast<uint8_t>((value >> (8 * i)) & 0xFF);
        }
    }

    /// Specialization for float (use memcpy to preserve IEEE 754 format)
    template <> inline void write_le<float>(std::vector<uint8_t> &buf, size_t &pos, float value) {
        uint32_t bits;
        std::memcpy(&bits, &value, sizeof(value));
        write_le(buf, pos, bits);
    }

    /// Specialization for double (use memcpy to preserve IEEE 754 format)
    template <> inline void write_le<double>(std::vector<uint8_t> &buf, size_t &pos, double value) {
        uint64_t bits;
        std::memcpy(&bits, &value, sizeof(value));
        write_le(buf, pos, bits);
    }

    /// Write options for GeoTIFF output
    struct WriteOptions {
        // No compression options - only uncompressed TIFF supported
    };

    /// Write out all layers in rc as a chained-IFD GeoTIFF.
    /// Each IFD can have its own CRS/DATUM/HEADING/PixelScale and custom tags.
    /// Supports multiple data types: uint8, int8, uint16, int16, uint32, int32, float, double
    /// Only uncompressed TIFF is supported
    ///
    /// CRS Flavor Handling:
    /// - ENU flavor: Grid data is already in local space, datum provides reference
    /// - WGS flavor: Grid data represents WGS coordinates, datum provides reference
    /// The Grid object contains the appropriate coordinate system based on parsing
    inline std::vector<uint8_t> toTiffBytes(RasterCollection const &rc, const WriteOptions &options = WriteOptions{}) {
        size_t N = rc.layers.size();
        if (N == 0)
            throw std::runtime_error("toTiffBytes(): no layers");

        // --- 1) Flatten each layer's grid into its own chunky strip ---
        // Now handles multiple data types via std::visit
        std::vector<std::vector<uint8_t>> strips(N);
        std::vector<uint32_t> stripCounts(N), stripOffsets(N);
        std::vector<uint16_t> bitsPerSample(N);
        std::vector<SampleFormat> sampleFormats(N);
        std::vector<uint16_t> samplesPerPixel(N);
        std::vector<PhotometricInterpretation> photometricInterp(N);

        // No compression support - data is stored uncompressed

        for (size_t i = 0; i < N; ++i) {
            auto const &layer = rc.layers[i];

            // Get bits per sample, sample format, samples per pixel, and photometric from the grid variant
            bitsPerSample[i] = get_bits_per_sample(layer.grid);
            sampleFormats[i] = get_sample_format(layer.grid);
            samplesPerPixel[i] = get_samples_per_pixel(layer.grid);
            photometricInterp[i] = get_photometric_interpretation(layer.grid);

            // Use std::visit to handle all grid types
            std::visit(
                [&](const auto &g) {
                    using GridType = std::decay_t<decltype(g)>;
                    using T = grid_element_type_t<GridType>;

                    uint32_t W = static_cast<uint32_t>(g.cols);
                    uint32_t H = static_cast<uint32_t>(g.rows);

                    std::vector<uint8_t> rawData;

                    if constexpr (std::is_same_v<T, RGBA>) {
                        // RGBA grid: write 4 bytes per pixel (R, G, B, A)
                        size_t sz = size_t(W) * H * 4;
                        rawData.resize(sz);

                        size_t idx = 0;
                        for (int32_t r = H - 1; r >= 0; --r) {
                            for (uint32_t c = 0; c < W; ++c) {
                                const RGBA &v = g(r, c);
                                rawData[idx++] = v.r;
                                rawData[idx++] = v.g;
                                rawData[idx++] = v.b;
                                rawData[idx++] = v.a;
                            }
                        }
                    } else {
                        // Scalar grid: write single value per pixel
                        uint32_t S = layer.samplesPerPixel > 0 ? layer.samplesPerPixel : 1;
                        size_t bytesPerSample = sizeof(T);
                        size_t sz = size_t(W) * H * S * bytesPerSample;

                        rawData.resize(sz);

                        // Fill chunky: band0,band1,... per pixel
                        // Flip vertically: write bottom row first to match GeoTIFF coordinate system
                        size_t idx = 0;
                        for (int32_t r = H - 1; r >= 0; --r) {
                            for (uint32_t c = 0; c < W; ++c) {
                                T v = g(r, c);
                                for (uint32_t s = 0; s < S; ++s) {
                                    write_le<T>(rawData, idx, v);
                                }
                            }
                        }
                    }

                    // Store uncompressed data directly
                    strips[i] = std::move(rawData);

                    stripCounts[i] = static_cast<uint32_t>(strips[i].size());
                },
                layer.grid);
        }

        // --- 2) Compute strip offsets (right after the 8-byte TIFF header) ---
        uint32_t p = 8;
        for (size_t i = 0; i < N; ++i) {
            stripOffsets[i] = p;
            p += stripCounts[i];
        }
        uint32_t firstIFD = p;

        // --- 3) Prepare per-layer metadata ---
        std::vector<std::string> descriptions(N);
        std::vector<uint32_t> descLengths(N);
        std::vector<uint32_t> descOffsets(N);
        std::vector<uint32_t> scaleOffsets(N);
        std::vector<uint32_t> geoKeyOffsets(N);
        std::vector<uint32_t> tiepointOffsets(N);
        std::vector<uint32_t> transformOffsets(N); // For ModelTransformationTag (rotated grids)
        std::vector<bool> layerHasRotation(N);     // Track which layers need transformation matrix

        for (size_t i = 0; i < N; ++i) {
            auto const &layer = rc.layers[i];
            layerHasRotation[i] = has_rotation(layer);

            // Build ImageDescription for this layer
            if (!layer.imageDescription.empty()) {
                descriptions[i] = layer.imageDescription; // Use custom description if provided
            } else {
                // Generate geospatial description (always WGS84)
                descriptions[i] = "CRS WGS84 DATUM " + std::to_string(layer.datum.latitude) + " " +
                                  std::to_string(layer.datum.longitude) + " " + std::to_string(layer.datum.altitude) +
                                  " SHIFT " + std::to_string(layer.shift.point.x) + " " +
                                  std::to_string(layer.shift.point.y) + " " + std::to_string(layer.shift.point.z) +
                                  " " + std::to_string(layer.shift.rotation.to_euler().yaw);
            }

            descLengths[i] = uint32_t(descriptions[i].size() + 1);
        }

        // --- 4) Compute IFD offsets and sizes ---
        std::vector<uint16_t> entryCounts(N);
        std::vector<uint32_t> customDataOffsets(N);
        std::vector<uint32_t> customDataSizes(N);

        for (size_t i = 0; i < N; ++i) {
            // Base tags: 9 standard + ImageDescription + PlanarConfig + SampleFormat + GeoKeyDirectory + custom tags
            // Plus either (ModelPixelScale + ModelTiepoint) OR ModelTransformation
            // Plus ExtraSamples tag if RGBA (4 samples per pixel)
            uint16_t geoTagCount = layerHasRotation[i] ? 1 : 2;          // 1 for transform, 2 for scale+tiepoint
            uint16_t extraSamplesTag = (samplesPerPixel[i] > 3) ? 1 : 0; // ExtraSamples for alpha channel
            entryCounts[i] = 9 + 1 + 1 + 1 + 1 + geoTagCount + extraSamplesTag +
                             static_cast<uint16_t>(rc.layers[i].customTags.size());

            // Calculate space needed for multi-value custom tag data
            customDataSizes[i] = 0;
            for (const auto &[tag, values] : rc.layers[i].customTags) {
                if (values.size() > 1) {
                    customDataSizes[i] += static_cast<uint32_t>(values.size() * 4); // 4 bytes per uint32_t
                }
            }
        }

        std::vector<uint32_t> ifdSizes(N);
        for (size_t i = 0; i < N; ++i) {
            ifdSizes[i] = 2 + entryCounts[i] * 12 + 4; // entry count + entries + next IFD pointer
        }

        std::vector<uint32_t> ifdOffsets(N);
        p = firstIFD;
        for (size_t i = 0; i < N; ++i) {
            ifdOffsets[i] = p;
            p += ifdSizes[i];
        }

        // --- 5) Compute offsets for variable-length data ---
        for (size_t i = 0; i < N; ++i) {
            descOffsets[i] = p;
            p += descLengths[i];

            geoKeyOffsets[i] = p;
            p += 40; // GeoKeyDirectory = 40 bytes (header:8 + 4 keys×8 = 40)

            if (layerHasRotation[i]) {
                // ModelTransformationTag: 16 doubles = 128 bytes
                transformOffsets[i] = p;
                p += 128;
                scaleOffsets[i] = 0;    // Not used
                tiepointOffsets[i] = 0; // Not used
            } else {
                // ModelPixelScaleTag + ModelTiepointTag
                scaleOffsets[i] = p;
                p += 24; // 3 doubles = 24 bytes

                tiepointOffsets[i] = p;
                p += 48; // 6 doubles = 48 bytes

                transformOffsets[i] = 0; // Not used
            }

            // Custom tag data offset
            customDataOffsets[i] = p;
            p += customDataSizes[i];
        }

        uint32_t totalSize = p;

        // --- 5) Allocate final buffer ---
        std::vector<uint8_t> buf(totalSize);
        size_t writePos = 0;

        // LE writers
        auto writeLE16 = [&](uint16_t v) {
            buf[writePos++] = uint8_t(v & 0xFF);
            buf[writePos++] = uint8_t(v >> 8);
        };
        auto writeLE32 = [&](uint32_t v) {
            buf[writePos++] = uint8_t(v & 0xFF);
            buf[writePos++] = uint8_t((v >> 8) & 0xFF);
            buf[writePos++] = uint8_t((v >> 16) & 0xFF);
            buf[writePos++] = uint8_t((v >> 24) & 0xFF);
        };
        auto writeDouble = [&](double d) {
            uint64_t bits;
            std::memcpy(&bits, &d, sizeof(d));
            for (int i = 0; i < 8; ++i) {
                buf[writePos++] = uint8_t((bits >> (8 * i)) & 0xFF);
            }
        };

        // --- 6) TIFF header ---
        buf[writePos++] = 'I';
        buf[writePos++] = 'I'; // little-endian
        writeLE16(42);         // magic
        writeLE32(firstIFD);   // offset to first IFD

        // --- 7) Pixel data strips ---
        for (size_t i = 0; i < N; ++i) {
            std::memcpy(&buf[writePos], strips[i].data(), strips[i].size());
            writePos += strips[i].size();
        }

        // --- 8) IFDs (with proper tag ordering per TIFF 6.0 spec) ---
        for (size_t i = 0; i < N; ++i) {
            // Seek to the correct IFD position
            writePos = ifdOffsets[i];

            auto const &layer = rc.layers[i];
            auto [rows, cols] = get_grid_dimensions(layer.grid);
            uint32_t W = static_cast<uint32_t>(cols);
            uint32_t H = static_cast<uint32_t>(rows);
            uint32_t S = samplesPerPixel[i]; // Use computed samples per pixel (4 for RGBA, 1 for scalar)
            uint32_t PC = layer.planarConfig > 0 ? layer.planarConfig : 1;
            uint32_t PI = static_cast<uint32_t>(photometricInterp[i]); // PhotometricInterpretation

            // Collect all tags into a vector for sorting
            std::vector<IfdEntry> entries;
            entries.reserve(entryCounts[i]);

            // Standard TIFF tags
            entries.push_back({256, 4, 1, W});                           // ImageWidth
            entries.push_back({257, 4, 1, H});                           // ImageLength
            entries.push_back({258, 3, 1, bitsPerSample[i]});            // BitsPerSample
            entries.push_back({259, 3, 1, 1});                           // Compression (1 = None)
            entries.push_back({262, 3, 1, PI});                          // PhotometricInterpretation
            entries.push_back({270, 2, descLengths[i], descOffsets[i]}); // ImageDescription
            entries.push_back({273, 4, 1, stripOffsets[i]});             // StripOffsets
            entries.push_back({277, 3, 1, S});                           // SamplesPerPixel
            entries.push_back({278, 4, 1, H});                           // RowsPerStrip
            entries.push_back({279, 4, 1, stripCounts[i]});              // StripByteCounts
            entries.push_back({284, 3, 1, PC});                          // PlanarConfiguration

            // ExtraSamples tag (338) for RGBA - indicates alpha channel type
            // Value 2 = Unassociated alpha (not pre-multiplied)
            if (S > 3) {
                entries.push_back({338, 3, 1, 2}); // ExtraSamples: unassociated alpha
            }

            entries.push_back({339, 3, 1, static_cast<uint32_t>(sampleFormats[i])}); // SampleFormat

            // GeoTIFF tags - use either Scale+Tiepoint or Transformation based on rotation
            if (layerHasRotation[i]) {
                // Rotated grid: use ModelTransformationTag (34264)
                entries.push_back({34264, 12, 16, transformOffsets[i]}); // ModelTransformationTag
            } else {
                // Non-rotated grid: use Scale + Tiepoint
                entries.push_back({33550, 12, 3, scaleOffsets[i]});    // ModelPixelScaleTag
                entries.push_back({33922, 12, 6, tiepointOffsets[i]}); // ModelTiepointTag
            }
            entries.push_back({34735, 3, 20, geoKeyOffsets[i]}); // GeoKeyDirectoryTag

            // Custom tags
            uint32_t customDataPos = customDataOffsets[i];
            for (const auto &[tag, values] : layer.customTags) {
                uint32_t valueOrOffset;
                if (values.size() == 1) {
                    valueOrOffset = values[0]; // Value fits in offset field
                } else {
                    valueOrOffset = customDataPos; // Pointer to data
                    customDataPos += static_cast<uint32_t>(values.size() * 4);
                }
                entries.push_back({tag, 4, static_cast<uint32_t>(values.size()), valueOrOffset});
            }

            // Sort entries by tag number (TIFF 6.0 requirement)
            std::sort(entries.begin(), entries.end());

            // Write entry count
            writeLE16(static_cast<uint16_t>(entries.size()));

            // Write sorted entries
            for (const auto &entry : entries) {
                writeLE16(entry.tag);
                writeLE16(entry.type);
                writeLE32(entry.count);
                writeLE32(entry.value_or_offset);
            }

            // next IFD pointer
            uint32_t next = (i + 1 < N ? ifdOffsets[i + 1] : 0);
            writeLE32(next);
        }

        // --- 9) Write variable-length data for each layer ---
        for (size_t i = 0; i < N; ++i) {
            auto const &layer = rc.layers[i];
            auto [gridRows, gridCols] = get_grid_dimensions(layer.grid);
            uint32_t W = static_cast<uint32_t>(gridCols);
            uint32_t H = static_cast<uint32_t>(gridRows);

            // Description text + NUL
            writePos = descOffsets[i];
            std::memcpy(&buf[writePos], descriptions[i].data(), descriptions[i].size());
            buf[writePos + descriptions[i].size()] = '\0';

            // GeoKeyDirectory for this layer (always WGS84)
            writePos = geoKeyOffsets[i];
            writeLE16(1); // KeyDirectoryVersion
            writeLE16(1); // KeyRevision
            writeLE16(0); // MinorRevision
            writeLE16(4); // NumberOfKeys

            // Key entry 1: GTModelTypeGeoKey
            writeLE16(1024); // GTModelTypeGeoKey
            writeLE16(0);    // TIFFTagLocation (0 means value is in ValueOffset)
            writeLE16(1);    // Count
            writeLE16(2);    // 2=Geographic (WGS84)

            // Key entry 2: GTRasterTypeGeoKey
            writeLE16(1025); // GTRasterTypeGeoKey
            writeLE16(0);    // TIFFTagLocation
            writeLE16(1);    // Count
            writeLE16(1);    // RasterPixelIsArea

            // Key entry 3: GeographicTypeGeoKey - EPSG:4326 for WGS84
            writeLE16(2048); // GeographicTypeGeoKey
            writeLE16(0);    // TIFFTagLocation
            writeLE16(1);    // Count
            writeLE16(4326); // EPSG:4326 (WGS84)

            // Key entry 4: GeogAngularUnitsGeoKey - degrees for WGS84
            writeLE16(2054); // GeogAngularUnitsGeoKey
            writeLE16(0);    // TIFFTagLocation
            writeLE16(1);    // Count
            writeLE16(9102); // 9102=degree

            // Calculate common values needed for georeferencing
            cc::frame::ENU center_enu{layer.shift.point.x, layer.shift.point.y, layer.shift.point.z, layer.datum};
            cc::earth::WGS center_wgs = cc::frame::to_wgs(center_enu);

            // Calculate pixel scale in degrees
            cc::frame::ENU east_point_enu{layer.shift.point.x + layer.resolution, layer.shift.point.y,
                                          layer.shift.point.z, layer.datum};
            cc::frame::ENU north_point_enu{layer.shift.point.x, layer.shift.point.y + layer.resolution,
                                           layer.shift.point.z, layer.datum};
            cc::earth::WGS east_point_wgs = cc::frame::to_wgs(east_point_enu);
            cc::earth::WGS north_point_wgs = cc::frame::to_wgs(north_point_enu);

            double scale_x = east_point_wgs.longitude - center_wgs.longitude; // degrees per pixel (longitude)
            double scale_y = north_point_wgs.latitude - center_wgs.latitude;  // degrees per pixel (latitude)

            // Calculate grid extents
            double grid_width_meters = W * layer.resolution;
            double grid_height_meters = H * layer.resolution;

            cc::frame::ENU grid_west{layer.shift.point.x - grid_width_meters / 2.0, layer.shift.point.y,
                                     layer.shift.point.z, layer.datum};
            cc::frame::ENU grid_north{layer.shift.point.x, layer.shift.point.y + grid_height_meters / 2.0,
                                      layer.shift.point.z, layer.datum};
            cc::earth::WGS west_wgs = cc::frame::to_wgs(grid_west);
            cc::earth::WGS north_wgs = cc::frame::to_wgs(grid_north);

            double half_width_deg = center_wgs.longitude - west_wgs.longitude;
            double half_height_deg = north_wgs.latitude - center_wgs.latitude;

            double top_left_lon = center_wgs.longitude - half_width_deg;
            double top_left_lat = center_wgs.latitude + half_height_deg;

            if (layerHasRotation[i]) {
                // Write ModelTransformationTag (34264) - 4x4 affine transformation matrix
                // Matrix format (row-major): [a b 0 d; e f 0 h; 0 0 1 0; 0 0 0 1]
                // Where: a,b,e,f encode rotation+scale, d,h encode translation
                //
                // For GeoTIFF, the transformation is:
                //   X_world = a * col + b * row + d
                //   Y_world = e * col + f * row + h
                //
                // With rotation angle θ (yaw):
                //   a = scale_x * cos(θ)
                //   b = -scale_y * sin(θ)  (note: scale_y is positive, Y increases northward)
                //   e = scale_x * sin(θ)
                //   f = scale_y * cos(θ)
                //   d = top_left_lon (adjusted for rotation)
                //   h = top_left_lat (adjusted for rotation)

                writePos = transformOffsets[i];
                double yaw = layer.shift.rotation.to_euler().yaw;
                double cos_yaw = std::cos(yaw);
                double sin_yaw = std::sin(yaw);

                // For rotated grids, we need to compute the top-left corner after rotation
                // The grid center stays the same, but corners rotate around it
                // Top-left corner in unrotated grid is at (-W/2, +H/2) pixels from center
                // After rotation: x' = x*cos - y*sin, y' = x*sin + y*cos
                double half_w_pix = W / 2.0;
                double half_h_pix = H / 2.0;

                // Rotated top-left corner offset in pixels
                double tl_col_offset = -half_w_pix * cos_yaw - half_h_pix * sin_yaw;
                double tl_row_offset = -half_w_pix * sin_yaw + half_h_pix * cos_yaw;

                // Convert to degrees and add to center
                double rotated_tl_lon = center_wgs.longitude + tl_col_offset * scale_x;
                double rotated_tl_lat = center_wgs.latitude + tl_row_offset * scale_y;

                // Row 1: [a, b, 0, d]
                writeDouble(scale_x * cos_yaw);  // a: X scale with rotation
                writeDouble(-scale_y * sin_yaw); // b: Y contribution to X (rotation)
                writeDouble(0.0);                // 0
                writeDouble(rotated_tl_lon);     // d: X translation (top-left longitude)

                // Row 2: [e, f, 0, h]
                writeDouble(scale_x * sin_yaw); // e: X contribution to Y (rotation)
                writeDouble(scale_y * cos_yaw); // f: Y scale with rotation (negative for south-down)
                writeDouble(0.0);               // 0
                writeDouble(rotated_tl_lat);    // h: Y translation (top-left latitude)

                // Row 3: [0, 0, 1, 0]
                writeDouble(0.0);
                writeDouble(0.0);
                writeDouble(1.0);
                writeDouble(0.0);

                // Row 4: [0, 0, 0, 1]
                writeDouble(0.0);
                writeDouble(0.0);
                writeDouble(0.0);
                writeDouble(1.0);
            } else {
                // Write ModelPixelScaleTag (33550) - 3 doubles: X, Y, Z scale
                writePos = scaleOffsets[i];
                writeDouble(scale_x);  // X scale in degrees (longitude, positive = eastward)
                writeDouble(-scale_y); // Y scale in degrees (negative = rows increase southward)
                writeDouble(0.0);      // Z scale

                // Write ModelTiepointTag (33922) - 6 doubles: I, J, K, X, Y, Z
                writePos = tiepointOffsets[i];
                writeDouble(0.0);                 // I: pixel column 0 (left edge)
                writeDouble(0.0);                 // J: pixel row 0 (top edge)
                writeDouble(0.0);                 // K: always 0 for 2D
                writeDouble(top_left_lon);        // X: longitude of top-left corner
                writeDouble(top_left_lat);        // Y: latitude of top-left corner
                writeDouble(center_wgs.altitude); // Z: altitude
            }

            // Write custom tag data for this layer
            writePos = customDataOffsets[i];
            for (const auto &[tag, values] : layer.customTags) {
                if (values.size() > 1) {
                    for (uint32_t value : values) {
                        writeLE32(value);
                    }
                }
            }
        }

        return buf;
    }

    /// Write a multi-IFD GeoTIFF to disk
    /// @param rc RasterCollection to write
    /// @param outPath Output file path
    /// @param options Write options (compression, predictor, etc.)
    inline void WriteRasterCollection(RasterCollection const &rc, fs::path const &outPath,
                                      const WriteOptions &options = WriteOptions{}) {
        auto bytes = toTiffBytes(rc, options);
        std::ofstream ofs(outPath, std::ios::binary);
        if (!ofs)
            throw std::runtime_error("cannot open " + outPath.string());
        ofs.write(reinterpret_cast<char *>(bytes.data()), bytes.size());
    }

    /// Write a datapod Layer as a GeoTIFF with multiple IFDs (one per Z layer)
    template <typename T>
    void WriteLayerCollection(const dp::Layer<T> &layer3d, const fs::path &outPath,
                              const dp::Geo &datum = dp::Geo{0.001, 0.001, 1.0}) {

        RasterCollection rc;
        rc.datum = datum;
        rc.shift = layer3d.pose;
        rc.resolution = layer3d.resolution;

        // Extract each Z layer as a separate IFD
        for (size_t layerIdx = 0; layerIdx < layer3d.layers; ++layerIdx) {
            // Extract 2D grid from 3D layer
            auto grid2d = layer3d.extract_grid(layerIdx);

            // Get Z coordinate for this layer
            auto centerPoint = layer3d.get_point(0, 0, layerIdx);
            double zCoord = centerPoint.z;

            // Convert to uint8_t if needed, preserve original shift but update Z coordinate
            dp::Pose layerPose = layer3d.pose;
            layerPose.point.z = zCoord; // Update Z to this layer's Z coordinate
            auto uint8Grid =
                dp::make_grid<uint8_t>(grid2d.rows, grid2d.cols, grid2d.resolution, true, layerPose, uint8_t{0});

            // Copy data with type conversion
            for (size_t r = 0; r < grid2d.rows; ++r) {
                for (size_t c = 0; c < grid2d.cols; ++c) {
                    uint8Grid(r, c) =
                        static_cast<uint8_t>(std::min(255.0, std::max(0.0, static_cast<double>(grid2d(r, c)))));
                }
            }

            // Create layer
            Layer layer;
            layer.grid = std::move(uint8Grid);
            layer.width = static_cast<uint32_t>(layer3d.cols);
            layer.height = static_cast<uint32_t>(layer3d.rows);
            layer.samplesPerPixel = 1;
            layer.planarConfig = 1;
            layer.datum = datum;
            layer.shift = layer3d.pose;            // Use the original layer's pose
            layer.resolution = layer3d.resolution; // Use the original layer's resolution

            // Set layer description with both geospatial data AND layer-specific info
            layer.imageDescription =
                "CRS WGS84 DATUM " + std::to_string(datum.latitude) + " " + std::to_string(datum.longitude) + " " +
                std::to_string(datum.altitude) + " SHIFT " + std::to_string(layer3d.pose.point.x) + " " +
                std::to_string(layer3d.pose.point.y) + " " + std::to_string(layer3d.pose.point.z) + " " +
                std::to_string(layer3d.pose.rotation.to_euler().yaw) + " " +
                "LayerHeight=" + std::to_string(layer3d.layer_height) + " " +
                "Resolution=" + std::to_string(layer3d.resolution) + " " +
                "LayerCount=" + std::to_string(layer3d.layers) + " " + "LayerIndex=" + std::to_string(layerIdx) + " " +
                "LayerZ=" + std::to_string(zCoord);

            rc.layers.push_back(std::move(layer));
        }

        // Write to file
        WriteRasterCollection(rc, outPath);
    }

} // namespace geotiv
