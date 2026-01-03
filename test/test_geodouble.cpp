#include <datapod/datapod.hpp>
namespace dp = datapod;
#include "geotiv/geotiv.hpp"
#include <cmath>
#include <doctest/doctest.h>
#include <filesystem>
#include <fstream>

TEST_CASE("GeoDoubleParamsTag (34736) support") {
    size_t rows = 10, cols = 10;
    double cellSize = 1.0;
    dp::Geo datum{48.0, 11.0, 500.0};
    dp::Pose shift{dp::Point{0, 0, 0}, dp::Quaternion::from_euler(0, 0, 0)};

    SUBCASE("Basic GeoDoubleParams round-trip") {
        auto grid = dp::make_grid<float>(rows, cols, cellSize, true, shift, 0.0f);

        geotiv::RasterCollection rc;
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = cellSize;

        geotiv::Layer layer;
        layer.grid = std::move(grid);
        layer.width = static_cast<uint32_t>(cols);
        layer.height = static_cast<uint32_t>(rows);
        layer.samplesPerPixel = 1;
        layer.planarConfig = 1;
        layer.datum = datum;
        layer.shift = shift;
        layer.resolution = cellSize;
        layer.geoDoubleParams = {1.0, 2.0, 3.0};

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_geodouble_basic.tif";
        geotiv::WriteRasterCollection(rc, testFile);

        // Read back and verify
        auto rc2 = geotiv::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);
        REQUIRE(rc2.layers[0].geoDoubleParams.size() == 3);
        CHECK(rc2.layers[0].geoDoubleParams[0] == doctest::Approx(1.0));
        CHECK(rc2.layers[0].geoDoubleParams[1] == doctest::Approx(2.0));
        CHECK(rc2.layers[0].geoDoubleParams[2] == doctest::Approx(3.0));

        std::filesystem::remove(testFile);
    }

    SUBCASE("High precision values") {
        auto grid = dp::make_grid<double>(rows, cols, cellSize, true, shift, 0.0);

        geotiv::RasterCollection rc;
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = cellSize;

        geotiv::Layer layer;
        layer.grid = std::move(grid);
        layer.width = static_cast<uint32_t>(cols);
        layer.height = static_cast<uint32_t>(rows);
        layer.samplesPerPixel = 1;
        layer.planarConfig = 1;
        layer.datum = datum;
        layer.shift = shift;
        layer.resolution = cellSize;
        // High precision values
        layer.geoDoubleParams = {
            3.141592653589793,    // Pi
            2.718281828459045,    // e
            1.414213562373095,    // sqrt(2)
            0.000000001234567890, // Very small
            123456789.987654321   // Large with decimals
        };

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_geodouble_precision.tif";
        geotiv::WriteRasterCollection(rc, testFile);

        // Read back and verify precision
        auto rc2 = geotiv::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);
        REQUIRE(rc2.layers[0].geoDoubleParams.size() == 5);
        CHECK(rc2.layers[0].geoDoubleParams[0] == doctest::Approx(3.141592653589793).epsilon(1e-15));
        CHECK(rc2.layers[0].geoDoubleParams[1] == doctest::Approx(2.718281828459045).epsilon(1e-15));
        CHECK(rc2.layers[0].geoDoubleParams[2] == doctest::Approx(1.414213562373095).epsilon(1e-15));
        CHECK(rc2.layers[0].geoDoubleParams[3] == doctest::Approx(0.000000001234567890).epsilon(1e-15));
        CHECK(rc2.layers[0].geoDoubleParams[4] == doctest::Approx(123456789.987654321).epsilon(1e-15));

        std::filesystem::remove(testFile);
    }

    SUBCASE("Empty GeoDoubleParams (tag should not be written)") {
        auto grid = dp::make_grid<uint8_t>(rows, cols, cellSize, true, shift, uint8_t{128});

        geotiv::RasterCollection rc;
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = cellSize;

        geotiv::Layer layer;
        layer.grid = std::move(grid);
        layer.width = static_cast<uint32_t>(cols);
        layer.height = static_cast<uint32_t>(rows);
        layer.samplesPerPixel = 1;
        layer.planarConfig = 1;
        layer.datum = datum;
        layer.shift = shift;
        layer.resolution = cellSize;
        // geoDoubleParams is empty (default)

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_geodouble_empty.tif";
        geotiv::WriteRasterCollection(rc, testFile);

        // Read back and verify tag is not present
        auto rc2 = geotiv::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);
        CHECK(rc2.layers[0].geoDoubleParams.empty());

        std::filesystem::remove(testFile);
    }

    SUBCASE("Single double value") {
        auto grid = dp::make_grid<int16_t>(rows, cols, cellSize, true, shift, int16_t{0});

        geotiv::RasterCollection rc;
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = cellSize;

        geotiv::Layer layer;
        layer.grid = std::move(grid);
        layer.width = static_cast<uint32_t>(cols);
        layer.height = static_cast<uint32_t>(rows);
        layer.samplesPerPixel = 1;
        layer.planarConfig = 1;
        layer.datum = datum;
        layer.shift = shift;
        layer.resolution = cellSize;
        layer.geoDoubleParams = {42.123456789};

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_geodouble_single.tif";
        geotiv::WriteRasterCollection(rc, testFile);

        // Read back and verify
        auto rc2 = geotiv::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);
        REQUIRE(rc2.layers[0].geoDoubleParams.size() == 1);
        CHECK(rc2.layers[0].geoDoubleParams[0] == doctest::Approx(42.123456789));

        std::filesystem::remove(testFile);
    }

    SUBCASE("Many double values") {
        auto grid = dp::make_grid<float>(rows, cols, cellSize, true, shift, 0.0f);

        geotiv::RasterCollection rc;
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = cellSize;

        geotiv::Layer layer;
        layer.grid = std::move(grid);
        layer.width = static_cast<uint32_t>(cols);
        layer.height = static_cast<uint32_t>(rows);
        layer.samplesPerPixel = 1;
        layer.planarConfig = 1;
        layer.datum = datum;
        layer.shift = shift;
        layer.resolution = cellSize;

        // Create array of 100 double values
        std::vector<double> values;
        for (int i = 0; i < 100; ++i) {
            values.push_back(i * 0.123456789);
        }
        layer.geoDoubleParams = values;

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_geodouble_many.tif";
        geotiv::WriteRasterCollection(rc, testFile);

        // Read back and verify
        auto rc2 = geotiv::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);
        REQUIRE(rc2.layers[0].geoDoubleParams.size() == 100);
        for (int i = 0; i < 100; ++i) {
            CHECK(rc2.layers[0].geoDoubleParams[i] == doctest::Approx(i * 0.123456789));
        }

        std::filesystem::remove(testFile);
    }

    SUBCASE("Negative and zero values") {
        auto grid = dp::make_grid<float>(rows, cols, cellSize, true, shift, 0.0f);

        geotiv::RasterCollection rc;
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = cellSize;

        geotiv::Layer layer;
        layer.grid = std::move(grid);
        layer.width = static_cast<uint32_t>(cols);
        layer.height = static_cast<uint32_t>(rows);
        layer.samplesPerPixel = 1;
        layer.planarConfig = 1;
        layer.datum = datum;
        layer.shift = shift;
        layer.resolution = cellSize;
        layer.geoDoubleParams = {-123.456, 0.0, -0.0, 456.789, -999.999};

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_geodouble_negative.tif";
        geotiv::WriteRasterCollection(rc, testFile);

        // Read back and verify
        auto rc2 = geotiv::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);
        REQUIRE(rc2.layers[0].geoDoubleParams.size() == 5);
        CHECK(rc2.layers[0].geoDoubleParams[0] == doctest::Approx(-123.456));
        CHECK(rc2.layers[0].geoDoubleParams[1] == doctest::Approx(0.0));
        CHECK(rc2.layers[0].geoDoubleParams[2] == doctest::Approx(-0.0));
        CHECK(rc2.layers[0].geoDoubleParams[3] == doctest::Approx(456.789));
        CHECK(rc2.layers[0].geoDoubleParams[4] == doctest::Approx(-999.999));

        std::filesystem::remove(testFile);
    }

    SUBCASE("Verify GeoDoubleParamsTag (34736) is written") {
        auto grid = dp::make_grid<float>(rows, cols, cellSize, true, shift, 0.0f);

        geotiv::RasterCollection rc;
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = cellSize;

        geotiv::Layer layer;
        layer.grid = std::move(grid);
        layer.width = static_cast<uint32_t>(cols);
        layer.height = static_cast<uint32_t>(rows);
        layer.samplesPerPixel = 1;
        layer.planarConfig = 1;
        layer.datum = datum;
        layer.shift = shift;
        layer.resolution = cellSize;
        layer.geoDoubleParams = {1.0, 2.0};

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_geodouble_tag34736.tif";
        geotiv::WriteRasterCollection(rc, testFile);

        // Parse IFD to verify tag 34736 exists
        std::ifstream f(testFile, std::ios::binary);
        REQUIRE(f.is_open());

        f.seekg(4);
        uint32_t ifdOffset = 0;
        f.read(reinterpret_cast<char *>(&ifdOffset), 4);

        f.seekg(ifdOffset);
        uint16_t numEntries = 0;
        f.read(reinterpret_cast<char *>(&numEntries), 2);

        bool found34736 = false;
        for (uint16_t e = 0; e < numEntries; ++e) {
            uint16_t tag = 0;
            f.read(reinterpret_cast<char *>(&tag), 2);
            if (tag == 34736) {
                found34736 = true;
                break;
            }
            f.seekg(10, std::ios::cur);
        }

        CHECK(found34736); // GeoDoubleParamsTag

        f.close();
        std::filesystem::remove(testFile);
    }

    SUBCASE("Special double values (infinity, very large/small)") {
        auto grid = dp::make_grid<float>(rows, cols, cellSize, true, shift, 0.0f);

        geotiv::RasterCollection rc;
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = cellSize;

        geotiv::Layer layer;
        layer.grid = std::move(grid);
        layer.width = static_cast<uint32_t>(cols);
        layer.height = static_cast<uint32_t>(rows);
        layer.samplesPerPixel = 1;
        layer.planarConfig = 1;
        layer.datum = datum;
        layer.shift = shift;
        layer.resolution = cellSize;
        layer.geoDoubleParams = {1e308,  // Very large
                                 1e-308, // Very small
                                 -1e308, // Very large negative
                                 std::numeric_limits<double>::max(), std::numeric_limits<double>::min()};

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_geodouble_special.tif";
        geotiv::WriteRasterCollection(rc, testFile);

        // Read back and verify
        auto rc2 = geotiv::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);
        REQUIRE(rc2.layers[0].geoDoubleParams.size() == 5);
        CHECK(rc2.layers[0].geoDoubleParams[0] == doctest::Approx(1e308));
        CHECK(rc2.layers[0].geoDoubleParams[1] == doctest::Approx(1e-308));
        CHECK(rc2.layers[0].geoDoubleParams[2] == doctest::Approx(-1e308));
        CHECK(rc2.layers[0].geoDoubleParams[3] == doctest::Approx(std::numeric_limits<double>::max()));
        CHECK(rc2.layers[0].geoDoubleParams[4] == doctest::Approx(std::numeric_limits<double>::min()));

        std::filesystem::remove(testFile);
    }
}
