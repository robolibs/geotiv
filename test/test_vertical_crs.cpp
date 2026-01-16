#include <datapod/datapod.hpp>
namespace dp = datapod;
#include "geotiv/geotiv.hpp"
#include <doctest/doctest.h>
#include <filesystem>

TEST_CASE("Vertical CRS GeoKeys support") {
    size_t rows = 10, cols = 10;
    double cellSize = 1.0;
    dp::Geo datum{48.0, 11.0, 500.0};
    dp::Pose shift{dp::Point{0, 0, 0}, dp::Quaternion::from_euler(0, 0, 0)};

    SUBCASE("WGS84 ellipsoid height (default behavior)") {
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
        layer.verticalDatum = 5030; // WGS84 ellipsoid
        layer.verticalUnits = 9001; // meters

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_vertical_wgs84.tif";
        geotiv::WriteRasterCollection(rc, testFile);

        // Read back and verify
        auto rc2 = geotiv::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);
        CHECK(rc2.layers[0].verticalDatum.has_value());
        CHECK(rc2.layers[0].verticalDatum.value() == 5030);
        CHECK(rc2.layers[0].verticalUnits.has_value());
        CHECK(rc2.layers[0].verticalUnits.value() == 9001);

        std::filesystem::remove(testFile);
    }

    SUBCASE("EGM96 geoid") {
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
        layer.verticalDatum = 5103; // EGM96 geoid
        layer.verticalUnits = 9001; // meters
        layer.verticalCitation = "EGM96 geoid height";

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_vertical_egm96.tif";
        geotiv::WriteRasterCollection(rc, testFile);

        // Read back and verify
        auto rc2 = geotiv::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);
        CHECK(rc2.layers[0].verticalDatum.has_value());
        CHECK(rc2.layers[0].verticalDatum.value() == 5103);
        CHECK(rc2.layers[0].verticalUnits.has_value());
        CHECK(rc2.layers[0].verticalUnits.value() == 9001);
        CHECK(rc2.layers[0].verticalCitation == "EGM96 geoid height");

        std::filesystem::remove(testFile);
    }

    SUBCASE("EGM2008 geoid") {
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
        layer.verticalDatum = 5773; // EGM2008 geoid
        layer.verticalUnits = 9001; // meters
        layer.verticalCitation = "EGM2008 geoid height";

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_vertical_egm2008.tif";
        geotiv::WriteRasterCollection(rc, testFile);

        // Read back and verify
        auto rc2 = geotiv::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);
        CHECK(rc2.layers[0].verticalDatum.has_value());
        CHECK(rc2.layers[0].verticalDatum.value() == 5773);
        CHECK(rc2.layers[0].verticalUnits.has_value());
        CHECK(rc2.layers[0].verticalUnits.value() == 9001);
        CHECK(rc2.layers[0].verticalCitation == "EGM2008 geoid height");

        std::filesystem::remove(testFile);
    }

    SUBCASE("Vertical units in feet") {
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
        layer.verticalDatum = 5030; // WGS84 ellipsoid
        layer.verticalUnits = 9002; // feet

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_vertical_feet.tif";
        geotiv::WriteRasterCollection(rc, testFile);

        // Read back and verify
        auto rc2 = geotiv::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);
        CHECK(rc2.layers[0].verticalDatum.has_value());
        CHECK(rc2.layers[0].verticalDatum.value() == 5030);
        CHECK(rc2.layers[0].verticalUnits.has_value());
        CHECK(rc2.layers[0].verticalUnits.value() == 9002);

        std::filesystem::remove(testFile);
    }

    SUBCASE("No vertical CRS (optional fields not set)") {
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
        // No vertical CRS fields set

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_no_vertical.tif";
        geotiv::WriteRasterCollection(rc, testFile);

        // Read back and verify
        auto rc2 = geotiv::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);
        CHECK_FALSE(rc2.layers[0].verticalDatum.has_value());
        CHECK_FALSE(rc2.layers[0].verticalUnits.has_value());
        CHECK(rc2.layers[0].verticalCitation.empty());

        std::filesystem::remove(testFile);
    }

    SUBCASE("Vertical citation only") {
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
        layer.verticalCitation = "Custom vertical datum";

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_vertical_citation_only.tif";
        geotiv::WriteRasterCollection(rc, testFile);

        // Read back and verify
        auto rc2 = geotiv::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);
        CHECK(rc2.layers[0].verticalCitation == "Custom vertical datum");

        std::filesystem::remove(testFile);
    }

    SUBCASE("Combined with horizontal citation") {
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
        layer.geoAsciiParams = "WGS 84";
        layer.verticalDatum = 5103; // EGM96
        layer.verticalUnits = 9001; // meters
        layer.verticalCitation = "EGM96 geoid";

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_vertical_combined.tif";
        geotiv::WriteRasterCollection(rc, testFile);

        // Read back and verify
        auto rc2 = geotiv::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);
        CHECK(rc2.layers[0].geoAsciiParams == "WGS 84");
        CHECK(rc2.layers[0].verticalDatum.has_value());
        CHECK(rc2.layers[0].verticalDatum.value() == 5103);
        CHECK(rc2.layers[0].verticalUnits.has_value());
        CHECK(rc2.layers[0].verticalUnits.value() == 9001);
        CHECK(rc2.layers[0].verticalCitation == "EGM96 geoid");

        std::filesystem::remove(testFile);
    }
}
