#include <datapod/datapod.hpp>
namespace dp = datapod;
#include "rastkit/rastkit.hpp"
#include <doctest/doctest.h>
#include <filesystem>
#include <fstream>

TEST_CASE("GeoAsciiParamsTag (34737) support") {
    size_t rows = 10, cols = 10;
    double cellSize = 1.0;
    dp::Geo datum{48.0, 11.0, 500.0};
    dp::Pose shift{dp::Point{0, 0, 0}, dp::Quaternion::from_euler(0, 0, 0)};

    SUBCASE("Basic GeoAsciiParams round-trip") {
        auto grid = dp::make_grid<float>(rows, cols, cellSize, true, shift, 0.0f);

        rastkit::RasterCollection rc;
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = cellSize;

        rastkit::Layer layer;
        layer.grid = std::move(grid);
        layer.width = static_cast<uint32_t>(cols);
        layer.height = static_cast<uint32_t>(rows);
        layer.samplesPerPixel = 1;
        layer.planarConfig = 1;
        layer.datum = datum;
        layer.shift = shift;
        layer.resolution = cellSize;
        layer.geoAsciiParams = "WGS 84";

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_geoascii_basic.tif";
        rastkit::WriteRasterCollection(rc, testFile);

        // Read back and verify
        auto rc2 = rastkit::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);
        CHECK(rc2.layers[0].geoAsciiParams == "WGS 84");

        std::filesystem::remove(testFile);
    }

    SUBCASE("Multiple strings with pipe delimiter") {
        auto grid = dp::make_grid<uint16_t>(rows, cols, cellSize, true, shift, uint16_t{100});

        rastkit::RasterCollection rc;
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = cellSize;

        rastkit::Layer layer;
        layer.grid = std::move(grid);
        layer.width = static_cast<uint32_t>(cols);
        layer.height = static_cast<uint32_t>(rows);
        layer.samplesPerPixel = 1;
        layer.planarConfig = 1;
        layer.datum = datum;
        layer.shift = shift;
        layer.resolution = cellSize;
        // Multiple strings separated by pipe
        layer.geoAsciiParams = "WGS 84|Geographic Coordinate System|EPSG:4326";

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_geoascii_multi.tif";
        rastkit::WriteRasterCollection(rc, testFile);

        // Read back and verify
        auto rc2 = rastkit::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);
        CHECK(rc2.layers[0].geoAsciiParams == "WGS 84|Geographic Coordinate System|EPSG:4326");

        std::filesystem::remove(testFile);
    }

    SUBCASE("Empty GeoAsciiParams (tag should not be written)") {
        auto grid = dp::make_grid<uint8_t>(rows, cols, cellSize, true, shift, uint8_t{128});

        rastkit::RasterCollection rc;
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = cellSize;

        rastkit::Layer layer;
        layer.grid = std::move(grid);
        layer.width = static_cast<uint32_t>(cols);
        layer.height = static_cast<uint32_t>(rows);
        layer.samplesPerPixel = 1;
        layer.planarConfig = 1;
        layer.datum = datum;
        layer.shift = shift;
        layer.resolution = cellSize;
        // geoAsciiParams is empty (default)

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_geoascii_empty.tif";
        rastkit::WriteRasterCollection(rc, testFile);

        // Read back and verify tag is not present
        auto rc2 = rastkit::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);
        CHECK(rc2.layers[0].geoAsciiParams.empty());

        std::filesystem::remove(testFile);
    }

    SUBCASE("Long citation string") {
        auto grid = dp::make_grid<double>(rows, cols, cellSize, true, shift, 123.456);

        rastkit::RasterCollection rc;
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = cellSize;

        rastkit::Layer layer;
        layer.grid = std::move(grid);
        layer.width = static_cast<uint32_t>(cols);
        layer.height = static_cast<uint32_t>(rows);
        layer.samplesPerPixel = 1;
        layer.planarConfig = 1;
        layer.datum = datum;
        layer.shift = shift;
        layer.resolution = cellSize;
        layer.geoAsciiParams =
            "WGS 84 / World Geodetic System 1984 - Geographic Coordinate System - EPSG:4326 - Horizontal datum: WGS84";

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_geoascii_long.tif";
        rastkit::WriteRasterCollection(rc, testFile);

        // Read back and verify
        auto rc2 = rastkit::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);
        CHECK(
            rc2.layers[0].geoAsciiParams ==
            "WGS 84 / World Geodetic System 1984 - Geographic Coordinate System - EPSG:4326 - Horizontal datum: WGS84");

        std::filesystem::remove(testFile);
    }

    SUBCASE("Verify GeoAsciiParamsTag (34737) is written") {
        auto grid = dp::make_grid<float>(rows, cols, cellSize, true, shift, 0.0f);

        rastkit::RasterCollection rc;
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = cellSize;

        rastkit::Layer layer;
        layer.grid = std::move(grid);
        layer.width = static_cast<uint32_t>(cols);
        layer.height = static_cast<uint32_t>(rows);
        layer.samplesPerPixel = 1;
        layer.planarConfig = 1;
        layer.datum = datum;
        layer.shift = shift;
        layer.resolution = cellSize;
        layer.geoAsciiParams = "Test Citation";

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_geoascii_tag34737.tif";
        rastkit::WriteRasterCollection(rc, testFile);

        // Parse IFD to verify tag 34737 exists
        std::ifstream f(testFile, std::ios::binary);
        REQUIRE(f.is_open());

        f.seekg(4);
        uint32_t ifdOffset = 0;
        f.read(reinterpret_cast<char *>(&ifdOffset), 4);

        f.seekg(ifdOffset);
        uint16_t numEntries = 0;
        f.read(reinterpret_cast<char *>(&numEntries), 2);

        bool found34737 = false;
        for (uint16_t e = 0; e < numEntries; ++e) {
            uint16_t tag = 0;
            f.read(reinterpret_cast<char *>(&tag), 2);
            if (tag == 34737) {
                found34737 = true;
                break;
            }
            f.seekg(10, std::ios::cur);
        }

        CHECK(found34737); // GeoAsciiParamsTag

        f.close();
        std::filesystem::remove(testFile);
    }

    SUBCASE("Special characters in GeoAsciiParams") {
        auto grid = dp::make_grid<int16_t>(rows, cols, cellSize, true, shift, int16_t{0});

        rastkit::RasterCollection rc;
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = cellSize;

        rastkit::Layer layer;
        layer.grid = std::move(grid);
        layer.width = static_cast<uint32_t>(cols);
        layer.height = static_cast<uint32_t>(rows);
        layer.samplesPerPixel = 1;
        layer.planarConfig = 1;
        layer.datum = datum;
        layer.shift = shift;
        layer.resolution = cellSize;
        // Test with special characters
        layer.geoAsciiParams = "WGS-84 (EPSG:4326) / Lat/Lon @ 0.5° resolution";

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_geoascii_special.tif";
        rastkit::WriteRasterCollection(rc, testFile);

        // Read back and verify
        auto rc2 = rastkit::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);
        CHECK(rc2.layers[0].geoAsciiParams == "WGS-84 (EPSG:4326) / Lat/Lon @ 0.5° resolution");

        std::filesystem::remove(testFile);
    }
}
