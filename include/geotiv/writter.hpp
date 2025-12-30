#pragma once

#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "concord/concord.hpp"
#include "geotiv/types.hpp"

namespace geotiv {

    namespace fs = std::filesystem;

    /// Write out all layers in rc as a chained-IFD GeoTIFF.
    /// Each IFD can have its own CRS/DATUM/HEADING/PixelScale and custom tags.
    ///
    /// CRS Flavor Handling:
    /// - ENU flavor: Grid data is already in local space, datum provides reference
    /// - WGS flavor: Grid data represents WGS coordinates, datum provides reference
    /// The Grid object contains the appropriate coordinate system based on parsing
    inline std::vector<uint8_t> toTiffBytes(RasterCollection const &rc) {
        size_t N = rc.layers.size();
        if (N == 0)
            throw std::runtime_error("toTiffBytes(): no layers");

        // --- 1) Flatten each layer's grid into its own chunky strip ---
        std::vector<std::vector<uint8_t>> strips(N);
        std::vector<uint32_t> stripCounts(N), stripOffsets(N);
        for (size_t i = 0; i < N; ++i) {
            auto const &layer = rc.layers[i];
            auto const &g = layer.grid;
            uint32_t W = static_cast<uint32_t>(g.cols());
            uint32_t H = static_cast<uint32_t>(g.rows());
            uint32_t S = layer.samplesPerPixel;
            size_t sz = size_t(W) * H * S;
            strips[i].resize(sz);
            stripCounts[i] = uint32_t(sz);
            // fill chunky: band0,band1,... per pixel
            // Flip vertically: write bottom row first to match GeoTIFF coordinate system
            size_t idx = 0;
            for (int32_t r = H - 1; r >= 0; --r) { // Start from bottom row, go up
                for (uint32_t c = 0; c < W; ++c) {
                    uint8_t v = g(r, c);
                    for (uint32_t s = 0; s < S; ++s) {
                        strips[i][idx++] = v;
                    }
                }
            }
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

        for (size_t i = 0; i < N; ++i) {
            auto const &layer = rc.layers[i];

            // Build ImageDescription for this layer
            if (!layer.imageDescription.empty()) {
                descriptions[i] = layer.imageDescription; // Use custom description if provided
            } else {
                // Generate geospatial description (always WGS84)
                descriptions[i] = "CRS WGS84 DATUM " + std::to_string(layer.datum.lat) + " " +
                                  std::to_string(layer.datum.lon) + " " + std::to_string(layer.datum.alt) + " SHIFT " +
                                  std::to_string(layer.shift.point.x) + " " + std::to_string(layer.shift.point.y) +
                                  " " + std::to_string(layer.shift.point.z) + " " +
                                  std::to_string(layer.shift.angle.yaw);
            }

            descLengths[i] = uint32_t(descriptions[i].size() + 1);
        }

        // --- 4) Compute IFD offsets and sizes ---
        std::vector<uint16_t> entryCounts(N);
        std::vector<uint32_t> customDataOffsets(N);
        std::vector<uint32_t> customDataSizes(N);

        for (size_t i = 0; i < N; ++i) {
            // Base tags: 9 standard + ImageDescription + PlanarConfig + ModelPixelScale + GeoKeyDirectory +
            // ModelTiepointTag + custom tags
            entryCounts[i] = 9 + 1 + 1 + 1 + 1 + 1 + static_cast<uint16_t>(rc.layers[i].customTags.size());

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

            scaleOffsets[i] = p;
            p += 24; // 3 doubles = 24 bytes

            geoKeyOffsets[i] = p;
            p += 40; // GeoKeyDirectory = 40 bytes (header:8 + 4 keys×8 = 40)

            tiepointOffsets[i] = p;
            p += 48; // ModelTiepointTag = 48 bytes (6 doubles)

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

        // --- 8) IFDs ---
        for (size_t i = 0; i < N; ++i) {
            // Seek to the correct IFD position
            writePos = ifdOffsets[i];

            writeLE16(entryCounts[i]);

            auto const &layer = rc.layers[i];
            auto const &g = layer.grid;
            uint32_t W = uint32_t(g.cols());
            uint32_t H = uint32_t(g.rows());
            uint32_t S = layer.samplesPerPixel;
            uint32_t PC = layer.planarConfig;

            // Tag 256: ImageWidth
            writeLE16(256);
            writeLE16(4);
            writeLE32(1);
            writeLE32(W);

            // Tag 257: ImageLength
            writeLE16(257);
            writeLE16(4);
            writeLE32(1);
            writeLE32(H);

            // Tag 258: BitsPerSample
            writeLE16(258);
            writeLE16(3);
            writeLE32(1);
            writeLE32(8);

            // Tag 259: Compression (1 = uncompressed)
            writeLE16(259);
            writeLE16(3);
            writeLE32(1);
            writeLE32(1);

            // Tag 262: PhotometricInterpretation (1 = BlackIsZero)
            writeLE16(262);
            writeLE16(3);
            writeLE32(1);
            writeLE32(1);

            // Tag 270: ImageDescription
            writeLE16(270);
            writeLE16(2);
            writeLE32(descLengths[i]);
            writeLE32(descOffsets[i]);

            // Tag 273: StripOffsets
            writeLE16(273);
            writeLE16(4);
            writeLE32(1);
            writeLE32(stripOffsets[i]);

            // Tag 277: SamplesPerPixel
            writeLE16(277);
            writeLE16(3);
            writeLE32(1);
            writeLE32(S);

            // Tag 278: RowsPerStrip
            writeLE16(278);
            writeLE16(4);
            writeLE32(1);
            writeLE32(H);

            // Tag 279: StripByteCounts
            writeLE16(279);
            writeLE16(4);
            writeLE32(1);
            writeLE32(stripCounts[i]);

            // Tag 284: PlanarConfiguration
            writeLE16(284);
            writeLE16(3);
            writeLE32(1);
            writeLE32(PC);

            // Tag 33550: ModelPixelScaleTag
            writeLE16(33550);
            writeLE16(12);
            writeLE32(3);
            writeLE32(scaleOffsets[i]);

            // Tag 33922: ModelTiepointTag (must come before 34735)
            writeLE16(33922);
            writeLE16(12);
            writeLE32(6);
            writeLE32(tiepointOffsets[i]);

            // Tag 34735: GeoKeyDirectoryTag
            writeLE16(34735);
            writeLE16(3);
            writeLE32(20); // 20 uint16 values: header(4) + 4 keys×4 = 20
            writeLE32(geoKeyOffsets[i]);

            // Write custom tags for this layer
            uint32_t customDataPos = customDataOffsets[i];
            for (const auto &[tag, values] : layer.customTags) {
                writeLE16(tag);
                writeLE16(4); // LONG type
                writeLE32(static_cast<uint32_t>(values.size()));
                if (values.size() == 1) {
                    writeLE32(values[0]); // Value fits in offset field
                } else {
                    writeLE32(customDataPos); // Pointer to data
                    customDataPos += static_cast<uint32_t>(values.size() * 4);
                }
            }

            // next IFD pointer
            uint32_t next = (i + 1 < N ? ifdOffsets[i + 1] : 0);
            writeLE32(next);
        }

        // --- 9) Write variable-length data for each layer ---
        for (size_t i = 0; i < N; ++i) {
            auto const &layer = rc.layers[i];

            // Description text + NUL
            writePos = descOffsets[i];
            std::memcpy(&buf[writePos], descriptions[i].data(), descriptions[i].size());
            buf[writePos + descriptions[i].size()] = '\0';

            // PixelScale doubles: X, Y, Z
            writePos = scaleOffsets[i];
            // Use precise concord library conversions for cm/mm accuracy
            // Create ENU points at grid center and at resolution distance from center
            concord::ENU center_enu{layer.shift.point.x, layer.shift.point.y, layer.shift.point.z, layer.datum};
            concord::ENU east_point_enu{layer.shift.point.x + layer.resolution, layer.shift.point.y,
                                        layer.shift.point.z, layer.datum}; // resolution meters east from center
            concord::ENU north_point_enu{layer.shift.point.x, layer.shift.point.y + layer.resolution,
                                         layer.shift.point.z, layer.datum}; // resolution meters north from center

            // Convert to WGS84 for precise degree differences
            concord::WGS center_wgs = center_enu.toWGS();
            concord::WGS east_point_wgs = east_point_enu.toWGS();
            concord::WGS north_point_wgs = north_point_enu.toWGS();

            // Calculate precise degree per meter scaling
            double resolution_deg_lon = east_point_wgs.lon - center_wgs.lon;  // precise longitude scale
            double resolution_deg_lat = north_point_wgs.lat - center_wgs.lat; // precise latitude scale

            writeDouble(resolution_deg_lon);  // X scale in degrees (longitude, positive = eastward)
            writeDouble(-resolution_deg_lat); // Y scale in degrees (negative = rows increase southward)
            writeDouble(0.0);                 // Z scale

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

            // Write ModelTiepointTag for this layer
            writePos = tiepointOffsets[i];
            auto const &g = layer.grid;
            uint32_t W = uint32_t(g.cols());
            uint32_t H = uint32_t(g.rows());

            // Tiepoint format: I,J,K,X,Y,Z where (I,J,K) are pixel coords and (X,Y,Z) are world coords
            // Standard practice: tie pixel (0,0) to top-left corner world coordinates
            writeDouble(0.0); // I: pixel column 0 (left edge)
            writeDouble(0.0); // J: pixel row 0 (top edge)
            writeDouble(0.0); // K: always 0 for 2D

            // Calculate top-left corner coordinates from grid center (shift)
            // First get grid dimensions in meters, then convert to degrees
            double grid_width_meters = W * layer.resolution;  // Total width in meters
            double grid_height_meters = H * layer.resolution; // Total height in meters

            // Convert half-extents from meters to degrees using precise conversion
            // We need to calculate how many degrees the grid spans from the grid center (shift point)
            concord::ENU grid_west{layer.shift.point.x - grid_width_meters / 2.0, layer.shift.point.y,
                                   layer.shift.point.z, layer.datum};
            concord::ENU grid_east{layer.shift.point.x + grid_width_meters / 2.0, layer.shift.point.y,
                                   layer.shift.point.z, layer.datum};
            concord::ENU grid_south{layer.shift.point.x, layer.shift.point.y - grid_height_meters / 2.0,
                                    layer.shift.point.z, layer.datum};
            concord::ENU grid_north{layer.shift.point.x, layer.shift.point.y + grid_height_meters / 2.0,
                                    layer.shift.point.z, layer.datum};

            concord::WGS west_wgs = grid_west.toWGS();
            concord::WGS east_wgs = grid_east.toWGS();
            concord::WGS south_wgs = grid_south.toWGS();
            concord::WGS north_wgs = grid_north.toWGS();

            double half_width = (east_wgs.lon - west_wgs.lon) / 2.0;    // Half width in degrees
            double half_height = (north_wgs.lat - south_wgs.lat) / 2.0; // Half height in degrees

            // Convert center to WGS84 first
            concord::ENU enuCenter{layer.shift.point.x, layer.shift.point.y, layer.shift.point.z, layer.datum};
            concord::WGS centerWGS = enuCenter.toWGS();

            // Calculate top-left corner from center
            double top_left_lon = centerWGS.lon - half_width;  // West of center
            double top_left_lat = centerWGS.lat + half_height; // North of center (top)

            writeDouble(top_left_lon);  // X: longitude of top-left corner
            writeDouble(top_left_lat);  // Y: latitude of top-left corner
            writeDouble(centerWGS.alt); // Z: altitude

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
    inline void WriteRasterCollection(RasterCollection const &rc, fs::path const &outPath) {
        auto bytes = toTiffBytes(rc);
        std::ofstream ofs(outPath, std::ios::binary);
        if (!ofs)
            throw std::runtime_error("cannot open " + outPath.string());
        ofs.write(reinterpret_cast<char *>(bytes.data()), bytes.size());
    }

    /// Write a Concord Layer as a GeoTIFF with multiple IFDs (one per Z layer)
    template <typename T>
    void WriteLayerCollection(const concord::Layer<T> &layer3d, const fs::path &outPath,
                              const concord::Datum &datum = concord::Datum{0.001, 0.001, 1.0}) {

        RasterCollection rc;
        rc.datum = datum;
        rc.shift = layer3d.shift();
        rc.resolution = layer3d.inradius();

        // Extract each Z layer as a separate IFD
        for (size_t layerIdx = 0; layerIdx < layer3d.layers(); ++layerIdx) {
            // Extract 2D grid from 3D layer
            auto grid2d = layer3d.extract_grid(layerIdx);

            // Get Z coordinate for this layer
            auto centerPoint = layer3d.get_point(0, 0, layerIdx);
            double zCoord = centerPoint.z;

            // Convert to uint8_t if needed, preserve original shift but update Z coordinate
            concord::Pose layerShift = layer3d.shift();
            layerShift.point.z = zCoord; // Update Z to this layer's Z coordinate
            concord::Grid<uint8_t> uint8Grid(grid2d.rows(), grid2d.cols(), grid2d.inradius(), true, layerShift);

            // Copy data with type conversion
            for (size_t r = 0; r < grid2d.rows(); ++r) {
                for (size_t c = 0; c < grid2d.cols(); ++c) {
                    uint8Grid(r, c) =
                        static_cast<uint8_t>(std::min(255.0, std::max(0.0, static_cast<double>(grid2d(r, c)))));
                }
            }

            // Create layer
            Layer layer;
            layer.grid = std::move(uint8Grid);
            layer.width = static_cast<uint32_t>(layer3d.cols());
            layer.height = static_cast<uint32_t>(layer3d.rows());
            layer.samplesPerPixel = 1;
            layer.planarConfig = 1;
            layer.datum = datum;
            layer.shift = layer3d.shift();         // Use the original layer's shift
            layer.resolution = layer3d.inradius(); // Use the original layer's resolution

            // Set layer description with both geospatial data AND layer-specific info
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

        // Write to file
        WriteRasterCollection(rc, outPath);
    }

} // namespace geotiv
