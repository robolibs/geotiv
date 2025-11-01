#pragma once

#include <filesystem>
#include <iosfwd>

#include "geotiv/types.hpp"

namespace geotiv {

    namespace detail {
        struct TIFFEntry {
            uint16_t tag, type;
            uint32_t count, valueOffset;
        };

        uint16_t readLE16(std::ifstream &f);
        uint32_t readLE32(std::ifstream &f);
        uint64_t readLE64(std::ifstream &f);
        uint16_t readBE16(std::ifstream &f);
        uint32_t readBE32(std::ifstream &f);
        uint64_t readBE64(std::ifstream &f);
        std::string readString(std::ifstream &f, uint32_t offset, uint32_t count);
    } // namespace detail

    RasterCollection ReadRasterCollection(const std::filesystem::path &file);

    std::ostream &operator<<(std::ostream &os, RasterCollection const &rc);

} // namespace geotiv
