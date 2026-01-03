#pragma once

#include <cstdint>
#include <sstream>
#include <stdexcept>
#include <string>

namespace geotiv {

    /// Base exception for all GeoTIFF errors
    class geotiff_error : public std::runtime_error {
      public:
        explicit geotiff_error(const std::string &message) : std::runtime_error(message) {}
    };

    /// Exception for file I/O errors
    class file_error : public geotiff_error {
      public:
        file_error(const std::string &filepath, const std::string &operation, const std::string &reason)
            : geotiff_error(format_message(filepath, operation, reason)), filepath_(filepath), operation_(operation),
              reason_(reason) {}

        const std::string &filepath() const { return filepath_; }
        const std::string &operation() const { return operation_; }
        const std::string &reason() const { return reason_; }

      private:
        std::string filepath_;
        std::string operation_;
        std::string reason_;

        static std::string format_message(const std::string &filepath, const std::string &operation,
                                          const std::string &reason) {
            std::ostringstream oss;
            oss << "File error in '" << filepath << "' during " << operation << ": " << reason;
            return oss.str();
        }
    };

    /// Exception for TIFF tag parsing errors
    class tag_error : public geotiff_error {
      public:
        tag_error(uint16_t tag, const std::string &tag_name, uint32_t offset, const std::string &reason)
            : geotiff_error(format_message(tag, tag_name, offset, reason)), tag_(tag), tag_name_(tag_name),
              offset_(offset), reason_(reason) {}

        uint16_t tag() const { return tag_; }
        const std::string &tag_name() const { return tag_name_; }
        uint32_t offset() const { return offset_; }
        const std::string &reason() const { return reason_; }

      private:
        uint16_t tag_;
        std::string tag_name_;
        uint32_t offset_;
        std::string reason_;

        static std::string format_message(uint16_t tag, const std::string &tag_name, uint32_t offset,
                                          const std::string &reason) {
            std::ostringstream oss;
            oss << "Tag " << tag << " (" << tag_name << ") at offset 0x" << std::hex << offset << std::dec << ": "
                << reason;
            return oss.str();
        }
    };

    /// Exception for validation errors
    class validation_error : public geotiff_error {
      public:
        validation_error(const std::string &field, const std::string &expected, const std::string &actual)
            : geotiff_error(format_message(field, expected, actual)), field_(field), expected_(expected),
              actual_(actual) {}

        const std::string &field() const { return field_; }
        const std::string &expected() const { return expected_; }
        const std::string &actual() const { return actual_; }

      private:
        std::string field_;
        std::string expected_;
        std::string actual_;

        static std::string format_message(const std::string &field, const std::string &expected,
                                          const std::string &actual) {
            std::ostringstream oss;
            oss << "Validation failed for " << field << ": expected " << expected << ", got " << actual;
            return oss.str();
        }
    };

    /// Exception for unsupported features
    class unsupported_error : public geotiff_error {
      public:
        unsupported_error(const std::string &feature, const std::string &supported_alternatives = "")
            : geotiff_error(format_message(feature, supported_alternatives)), feature_(feature),
              supported_alternatives_(supported_alternatives) {}

        const std::string &feature() const { return feature_; }
        const std::string &supported_alternatives() const { return supported_alternatives_; }

      private:
        std::string feature_;
        std::string supported_alternatives_;

        static std::string format_message(const std::string &feature, const std::string &supported_alternatives) {
            std::ostringstream oss;
            oss << "Unsupported feature: " << feature;
            if (!supported_alternatives.empty()) {
                oss << ". Supported: " << supported_alternatives;
            }
            return oss.str();
        }
    };

    /// Helper function to get tag name from tag number
    inline std::string get_tag_name(uint16_t tag) {
        switch (tag) {
        case 254:
            return "NewSubfileType";
        case 256:
            return "ImageWidth";
        case 257:
            return "ImageLength";
        case 258:
            return "BitsPerSample";
        case 259:
            return "Compression";
        case 262:
            return "PhotometricInterpretation";
        case 270:
            return "ImageDescription";
        case 273:
            return "StripOffsets";
        case 277:
            return "SamplesPerPixel";
        case 278:
            return "RowsPerStrip";
        case 279:
            return "StripByteCounts";
        case 282:
            return "XResolution";
        case 283:
            return "YResolution";
        case 284:
            return "PlanarConfiguration";
        case 296:
            return "ResolutionUnit";
        case 305:
            return "Software";
        case 306:
            return "DateTime";
        case 315:
            return "Artist";
        case 338:
            return "ExtraSamples";
        case 339:
            return "SampleFormat";
        case 33550:
            return "ModelPixelScaleTag";
        case 33922:
            return "ModelTiepointTag";
        case 34264:
            return "ModelTransformationTag";
        case 34735:
            return "GeoKeyDirectoryTag";
        case 34736:
            return "GeoDoubleParamsTag";
        case 34737:
            return "GeoAsciiParamsTag";
        case 42113:
            return "GDAL_NODATA";
        default:
            return "Unknown";
        }
    }

    /// Helper function to get compression name
    inline std::string get_compression_name(uint16_t compression) {
        switch (compression) {
        case 1:
            return "None (uncompressed)";
        case 2:
            return "CCITT 1D";
        case 3:
            return "Group 3 Fax";
        case 4:
            return "Group 4 Fax";
        case 5:
            return "LZW";
        case 6:
            return "JPEG (old-style)";
        case 7:
            return "JPEG";
        case 8:
            return "Deflate (Adobe-style)";
        case 32773:
            return "PackBits";
        case 32946:
            return "Deflate (PKZIP-style)";
        default:
            return "Unknown (" + std::to_string(compression) + ")";
        }
    }

} // namespace geotiv
