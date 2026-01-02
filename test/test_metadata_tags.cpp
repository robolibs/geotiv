#include <datapod/datapod.hpp>
namespace dp = datapod;
#include "geotiv/geotiv.hpp"
#include <doctest/doctest.h>
#include <filesystem>
#include <fstream>

TEST_CASE("Metadata tags (Software, DateTime, Resolution)") {
    size_t rows = 5, cols = 5;
    double cellSize = 1.0;
    dp::Geo datum{48.0, 11.0, 500.0};
    dp::Pose shift{dp::Point{0, 0, 0}, dp::Quaternion::from_euler(0, 0, 0)};

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

    rc.layers.push_back(std::move(layer));

    SUBCASE("Software tag (305) with default value") {
        auto bytes = geotiv::toTiffBytes(rc);
        std::string bytesStr(bytes.begin(), bytes.end());
        CHECK(bytesStr.find("geotiv 0.0.2") != std::string::npos);
    }

    SUBCASE("Software tag (305) with custom value") {
        geotiv::WriteOptions opts;
        opts.software = "My Custom GeoTIFF Writer v2.0";
        auto bytes = geotiv::toTiffBytes(rc, opts);
        std::string bytesStr(bytes.begin(), bytes.end());
        CHECK(bytesStr.find("My Custom GeoTIFF Writer v2.0") != std::string::npos);
    }

    SUBCASE("DateTime tag (306) with custom value") {
        geotiv::WriteOptions opts;
        opts.datetime = "2026:01:02 17:30:00";
        auto bytes = geotiv::toTiffBytes(rc, opts);
        std::string bytesStr(bytes.begin(), bytes.end());
        CHECK(bytesStr.find("2026:01:02 17:30:00") != std::string::npos);
    }

    SUBCASE("DateTime tag (306) with auto-generated value") {
        geotiv::WriteOptions opts;
        // datetime is empty, should auto-generate
        auto bytes = geotiv::toTiffBytes(rc, opts);
        std::string bytesStr(bytes.begin(), bytes.end());
        // Should contain year "202" (2020s)
        CHECK(bytesStr.find("202") != std::string::npos);
    }

    SUBCASE("Resolution tags (282/283/296) with default 72 DPI") {
        std::string testFile = "test_metadata_default_res.tif";
        geotiv::WriteRasterCollection(rc, testFile);

        // Parse IFD to verify tags exist
        std::ifstream f(testFile, std::ios::binary);
        REQUIRE(f.is_open());

        f.seekg(4);
        uint32_t ifdOffset = 0;
        f.read(reinterpret_cast<char *>(&ifdOffset), 4);

        f.seekg(ifdOffset);
        uint16_t numEntries = 0;
        f.read(reinterpret_cast<char *>(&numEntries), 2);

        bool found282 = false, found283 = false, found296 = false;
        for (uint16_t e = 0; e < numEntries; ++e) {
            uint16_t tag = 0;
            f.read(reinterpret_cast<char *>(&tag), 2);
            if (tag == 282)
                found282 = true;
            if (tag == 283)
                found283 = true;
            if (tag == 296)
                found296 = true;
            f.seekg(10, std::ios::cur);
        }

        CHECK(found282); // XResolution
        CHECK(found283); // YResolution
        CHECK(found296); // ResolutionUnit

        f.close();
        std::filesystem::remove(testFile);
    }

    SUBCASE("Resolution tags (282/283/296) with custom 300 DPI") {
        geotiv::WriteOptions opts;
        opts.xresolution_num = 300;
        opts.yresolution_num = 300;
        opts.resolution_unit = 2; // Inch

        std::string testFile = "test_metadata_300dpi.tif";
        geotiv::WriteRasterCollection(rc, testFile, opts);

        CHECK(std::filesystem::exists(testFile));
        std::filesystem::remove(testFile);
    }

    SUBCASE("All metadata tags together") {
        geotiv::WriteOptions opts;
        opts.software = "Test Suite v1.0";
        opts.datetime = "2026:01:02 18:00:00";
        opts.xresolution_num = 150;
        opts.yresolution_num = 150;
        opts.resolution_unit = 2;

        std::string testFile = "test_all_metadata.tif";
        geotiv::WriteRasterCollection(rc, testFile, opts);

        // Read back and verify
        auto rc2 = geotiv::ReadRasterCollection(testFile);
        CHECK(rc2.layers.size() == 1);

        // Verify data integrity
        const auto &grid2 = rc2.layers[0].gridAs<uint8_t>();
        CHECK(grid2.rows == rows);
        CHECK(grid2.cols == cols);

        std::filesystem::remove(testFile);
    }
}

TEST_CASE("NoData value support (GDAL_NODATA tag 42113)") {
    size_t rows = 10, cols = 10;
    double cellSize = 1.0;
    dp::Geo datum{48.0, 11.0, 500.0};
    dp::Pose shift{dp::Point{0, 0, 0}, dp::Quaternion::from_euler(0, 0, 0)};

    SUBCASE("NoData value with float32 data") {
        auto grid = dp::make_grid<float>(rows, cols, cellSize, true, shift, 0.0f);

        // Set some pixels to nodata value
        float nodataVal = -9999.0f;
        grid.at(0, 0) = nodataVal;
        grid.at(5, 5) = nodataVal;
        grid.at(9, 9) = nodataVal;

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
        layer.noDataValue = nodataVal; // Set nodata value

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_nodata_float.tif";
        geotiv::WriteRasterCollection(rc, testFile);

        // Read back and verify nodata value is preserved
        auto rc2 = geotiv::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);

        CHECK(rc2.layers[0].noDataValue.has_value());
        CHECK(rc2.layers[0].noDataValue.value() == doctest::Approx(nodataVal));

        // Verify data integrity
        const auto &grid2 = rc2.layers[0].gridAs<float>();
        CHECK(grid2.at(0, 0) == doctest::Approx(nodataVal));
        CHECK(grid2.at(5, 5) == doctest::Approx(nodataVal));
        CHECK(grid2.at(9, 9) == doctest::Approx(nodataVal));

        std::filesystem::remove(testFile);
    }

    SUBCASE("NoData value with int16 data") {
        auto grid = dp::make_grid<int16_t>(rows, cols, cellSize, true, shift, int16_t{100});

        int16_t nodataVal = -32768; // Common nodata for int16
        grid.at(0, 0) = nodataVal;
        grid.at(3, 7) = nodataVal;

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
        layer.noDataValue = static_cast<double>(nodataVal);

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_nodata_int16.tif";
        geotiv::WriteRasterCollection(rc, testFile);

        // Read back and verify
        auto rc2 = geotiv::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);

        CHECK(rc2.layers[0].noDataValue.has_value());
        CHECK(rc2.layers[0].noDataValue.value() == doctest::Approx(static_cast<double>(nodataVal)));

        const auto &grid2 = rc2.layers[0].gridAs<int16_t>();
        CHECK(grid2.at(0, 0) == nodataVal);
        CHECK(grid2.at(3, 7) == nodataVal);

        std::filesystem::remove(testFile);
    }

    SUBCASE("Layer without NoData value (tag should not be written)") {
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
        // noDataValue is NOT set (std::nullopt)

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_no_nodata.tif";
        geotiv::WriteRasterCollection(rc, testFile);

        // Read back and verify nodata is not present
        auto rc2 = geotiv::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);

        CHECK_FALSE(rc2.layers[0].noDataValue.has_value());

        std::filesystem::remove(testFile);
    }

    SUBCASE("NoData with double precision value") {
        auto grid = dp::make_grid<double>(rows, cols, cellSize, true, shift, 123.456);

        double nodataVal = -9999.123456789;
        grid.at(2, 2) = nodataVal;

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
        layer.noDataValue = nodataVal;

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_nodata_double.tif";
        geotiv::WriteRasterCollection(rc, testFile);

        // Read back and verify precision is maintained
        auto rc2 = geotiv::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);

        CHECK(rc2.layers[0].noDataValue.has_value());
        // Check with reasonable precision (string conversion may lose some precision)
        CHECK(rc2.layers[0].noDataValue.value() == doctest::Approx(nodataVal).epsilon(1e-9));

        std::filesystem::remove(testFile);
    }

    SUBCASE("Verify GDAL_NODATA tag (42113) is written") {
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
        layer.noDataValue = -9999.0;

        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_nodata_tag42113.tif";
        geotiv::WriteRasterCollection(rc, testFile);

        // Parse IFD to verify tag 42113 exists
        std::ifstream f(testFile, std::ios::binary);
        REQUIRE(f.is_open());

        f.seekg(4);
        uint32_t ifdOffset = 0;
        f.read(reinterpret_cast<char *>(&ifdOffset), 4);

        f.seekg(ifdOffset);
        uint16_t numEntries = 0;
        f.read(reinterpret_cast<char *>(&numEntries), 2);

        bool found42113 = false;
        for (uint16_t e = 0; e < numEntries; ++e) {
            uint16_t tag = 0;
            f.read(reinterpret_cast<char *>(&tag), 2);
            if (tag == 42113) {
                found42113 = true;
                break;
            }
            f.seekg(10, std::ios::cur);
        }

        CHECK(found42113); // GDAL_NODATA tag

        f.close();
        std::filesystem::remove(testFile);
    }
}
