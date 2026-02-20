#include <datapod/datapod.hpp>
namespace dp = datapod;
#include "rastkit/rastkit.hpp"
#include <cmath>
#include <doctest/doctest.h>
#include <filesystem>

TEST_CASE("GridVariant type support") {
    SUBCASE("uint8_t grid (default)") {
        auto grid = dp::make_grid<uint8_t>(10, 10, 1.0, true, dp::Pose{}, uint8_t{0});
        grid(5, 5) = 255;

        rastkit::GridVariant variant = std::move(grid);

        CHECK(rastkit::get_bits_per_sample(variant) == 8);
        CHECK(rastkit::get_sample_format(variant) == rastkit::SampleFormat::UnsignedInt);
        CHECK(rastkit::holds_grid_type<uint8_t>(variant));
        CHECK_FALSE(rastkit::holds_grid_type<uint16_t>(variant));

        auto &g = rastkit::get_grid<uint8_t>(variant);
        CHECK(g(5, 5) == 255);
    }

    SUBCASE("int8_t grid") {
        auto grid = dp::make_grid<int8_t>(10, 10, 1.0, true, dp::Pose{}, int8_t{0});
        grid(5, 5) = -128;
        grid(6, 6) = 127;

        rastkit::GridVariant variant = std::move(grid);

        CHECK(rastkit::get_bits_per_sample(variant) == 8);
        CHECK(rastkit::get_sample_format(variant) == rastkit::SampleFormat::SignedInt);
        CHECK(rastkit::holds_grid_type<int8_t>(variant));

        auto &g = rastkit::get_grid<int8_t>(variant);
        CHECK(g(5, 5) == -128);
        CHECK(g(6, 6) == 127);
    }

    SUBCASE("uint16_t grid") {
        auto grid = dp::make_grid<uint16_t>(10, 10, 1.0, true, dp::Pose{}, uint16_t{0});
        grid(5, 5) = 65535;
        grid(6, 6) = 1000;

        rastkit::GridVariant variant = std::move(grid);

        CHECK(rastkit::get_bits_per_sample(variant) == 16);
        CHECK(rastkit::get_sample_format(variant) == rastkit::SampleFormat::UnsignedInt);
        CHECK(rastkit::holds_grid_type<uint16_t>(variant));

        auto &g = rastkit::get_grid<uint16_t>(variant);
        CHECK(g(5, 5) == 65535);
        CHECK(g(6, 6) == 1000);
    }

    SUBCASE("int16_t grid") {
        auto grid = dp::make_grid<int16_t>(10, 10, 1.0, true, dp::Pose{}, int16_t{0});
        grid(5, 5) = -32768;
        grid(6, 6) = 32767;

        rastkit::GridVariant variant = std::move(grid);

        CHECK(rastkit::get_bits_per_sample(variant) == 16);
        CHECK(rastkit::get_sample_format(variant) == rastkit::SampleFormat::SignedInt);
        CHECK(rastkit::holds_grid_type<int16_t>(variant));

        auto &g = rastkit::get_grid<int16_t>(variant);
        CHECK(g(5, 5) == -32768);
        CHECK(g(6, 6) == 32767);
    }

    SUBCASE("uint32_t grid") {
        auto grid = dp::make_grid<uint32_t>(10, 10, 1.0, true, dp::Pose{}, uint32_t{0});
        grid(5, 5) = 4294967295U;
        grid(6, 6) = 1000000;

        rastkit::GridVariant variant = std::move(grid);

        CHECK(rastkit::get_bits_per_sample(variant) == 32);
        CHECK(rastkit::get_sample_format(variant) == rastkit::SampleFormat::UnsignedInt);
        CHECK(rastkit::holds_grid_type<uint32_t>(variant));

        auto &g = rastkit::get_grid<uint32_t>(variant);
        CHECK(g(5, 5) == 4294967295U);
        CHECK(g(6, 6) == 1000000);
    }

    SUBCASE("int32_t grid") {
        auto grid = dp::make_grid<int32_t>(10, 10, 1.0, true, dp::Pose{}, int32_t{0});
        grid(5, 5) = -2147483648;
        grid(6, 6) = 2147483647;

        rastkit::GridVariant variant = std::move(grid);

        CHECK(rastkit::get_bits_per_sample(variant) == 32);
        CHECK(rastkit::get_sample_format(variant) == rastkit::SampleFormat::SignedInt);
        CHECK(rastkit::holds_grid_type<int32_t>(variant));

        auto &g = rastkit::get_grid<int32_t>(variant);
        CHECK(g(5, 5) == -2147483648);
        CHECK(g(6, 6) == 2147483647);
    }

    SUBCASE("float grid") {
        auto grid = dp::make_grid<float>(10, 10, 1.0, true, dp::Pose{}, 0.0f);
        grid(5, 5) = 3.14159f;
        grid(6, 6) = -1.5e10f;

        rastkit::GridVariant variant = std::move(grid);

        CHECK(rastkit::get_bits_per_sample(variant) == 32);
        CHECK(rastkit::get_sample_format(variant) == rastkit::SampleFormat::Float);
        CHECK(rastkit::holds_grid_type<float>(variant));

        auto &g = rastkit::get_grid<float>(variant);
        CHECK(g(5, 5) == doctest::Approx(3.14159f));
        CHECK(g(6, 6) == doctest::Approx(-1.5e10f));
    }

    SUBCASE("double grid") {
        auto grid = dp::make_grid<double>(10, 10, 1.0, true, dp::Pose{}, 0.0);
        grid(5, 5) = 3.141592653589793;
        grid(6, 6) = -1.5e100;

        rastkit::GridVariant variant = std::move(grid);

        CHECK(rastkit::get_bits_per_sample(variant) == 64);
        CHECK(rastkit::get_sample_format(variant) == rastkit::SampleFormat::Float);
        CHECK(rastkit::holds_grid_type<double>(variant));

        auto &g = rastkit::get_grid<double>(variant);
        CHECK(g(5, 5) == doctest::Approx(3.141592653589793));
        CHECK(g(6, 6) == doctest::Approx(-1.5e100));
    }

    SUBCASE("get_grid_dimensions helper") {
        auto grid = dp::make_grid<uint16_t>(20, 30, 1.0, true, dp::Pose{}, uint16_t{0});
        rastkit::GridVariant variant = std::move(grid);

        auto [rows, cols] = rastkit::get_grid_dimensions(variant);
        CHECK(rows == 20);
        CHECK(cols == 30);
    }

    SUBCASE("get_grid_resolution helper") {
        auto grid = dp::make_grid<float>(10, 10, 2.5, true, dp::Pose{}, 0.0f);
        rastkit::GridVariant variant = std::move(grid);

        CHECK(rastkit::get_grid_resolution(variant) == doctest::Approx(2.5));
    }

    SUBCASE("RGBA grid") {
        auto grid = dp::make_grid<rastkit::RGBA>(10, 10, 1.0, true, dp::Pose{}, rastkit::RGBA{});
        grid(5, 5) = rastkit::RGBA{255, 128, 64, 200};
        grid(6, 6) = rastkit::RGBA{0, 0, 0, 255};

        rastkit::GridVariant variant = std::move(grid);

        CHECK(rastkit::get_bits_per_sample(variant) == 8);
        CHECK(rastkit::get_sample_format(variant) == rastkit::SampleFormat::UnsignedInt);
        CHECK(rastkit::get_samples_per_pixel(variant) == 4);
        CHECK(rastkit::get_photometric_interpretation(variant) == rastkit::PhotometricInterpretation::RGB);
        CHECK(rastkit::is_color_grid(variant));
        CHECK(rastkit::holds_grid_type<rastkit::RGBA>(variant));

        auto &g = rastkit::get_grid<rastkit::RGBA>(variant);
        CHECK(g(5, 5).r == 255);
        CHECK(g(5, 5).g == 128);
        CHECK(g(5, 5).b == 64);
        CHECK(g(5, 5).a == 200);
    }
}

TEST_CASE("Layer GridVariant convenience methods") {
    SUBCASE("Layer with uint8_t grid") {
        rastkit::Layer layer;
        auto grid = dp::make_grid<uint8_t>(10, 10, 1.0, true, dp::Pose{}, uint8_t{128});
        layer.grid = std::move(grid);

        CHECK(layer.bitsPerSample() == 8);
        CHECK(layer.sampleFormat() == rastkit::SampleFormat::UnsignedInt);
        CHECK(layer.holdsType<uint8_t>());
        CHECK_FALSE(layer.holdsType<float>());

        auto &g = layer.gridAs<uint8_t>();
        CHECK(g(0, 0) == 128);
    }

    SUBCASE("Layer with float grid") {
        rastkit::Layer layer;
        auto grid = dp::make_grid<float>(10, 10, 1.0, true, dp::Pose{}, 1.5f);
        layer.grid = std::move(grid);

        CHECK(layer.bitsPerSample() == 32);
        CHECK(layer.sampleFormat() == rastkit::SampleFormat::Float);
        CHECK(layer.holdsType<float>());

        auto *g = layer.gridIf<float>();
        REQUIRE(g != nullptr);
        CHECK((*g)(0, 0) == doctest::Approx(1.5f));

        auto *wrong = layer.gridIf<uint8_t>();
        CHECK(wrong == nullptr);
    }

    SUBCASE("Layer with RGBA grid") {
        rastkit::Layer layer;
        auto grid = dp::make_grid<rastkit::RGBA>(10, 10, 1.0, true, dp::Pose{}, rastkit::RGBA{255, 128, 64, 200});
        layer.grid = std::move(grid);

        CHECK(layer.bitsPerSample() == 8);
        CHECK(layer.sampleFormat() == rastkit::SampleFormat::UnsignedInt);
        CHECK(layer.photometricInterpretation() == rastkit::PhotometricInterpretation::RGB);
        CHECK(layer.isColorLayer());
        CHECK(layer.holdsType<rastkit::RGBA>());
        CHECK_FALSE(layer.holdsType<uint8_t>());

        auto *g = layer.gridIf<rastkit::RGBA>();
        REQUIRE(g != nullptr);
        CHECK((*g)(0, 0).r == 255);
        CHECK((*g)(0, 0).g == 128);
        CHECK((*g)(0, 0).b == 64);
        CHECK((*g)(0, 0).a == 200);
    }
}

TEST_CASE("Round-trip write/read for different grid types") {
    dp::Geo datum{48.0, 11.0, 500.0};
    dp::Pose shift{dp::Point{100, 200, 0}, dp::Quaternion::from_euler(0, 0, 0)};
    double resolution = 1.0;

    SUBCASE("uint16_t round-trip") {
        auto grid = dp::make_grid<uint16_t>(5, 5, resolution, true, shift, uint16_t{0});
        grid(0, 0) = 0;
        grid(2, 2) = 32768;
        grid(4, 4) = 65535;

        rastkit::RasterCollection rc;
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = resolution;

        rastkit::Layer layer;
        layer.grid = std::move(grid);
        layer.width = 5;
        layer.height = 5;
        layer.samplesPerPixel = 1;
        layer.planarConfig = 1;
        layer.datum = datum;
        layer.shift = shift;
        layer.resolution = resolution;
        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_uint16.tif";
        rastkit::WriteRasterCollection(rc, testFile);

        auto rc2 = rastkit::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);
        CHECK(rc2.layers[0].holdsType<uint16_t>());

        auto &g = rc2.layers[0].gridAs<uint16_t>();
        CHECK(g(0, 0) == 0);
        CHECK(g(2, 2) == 32768);
        CHECK(g(4, 4) == 65535);

        std::filesystem::remove(testFile);
    }

    SUBCASE("int16_t round-trip") {
        auto grid = dp::make_grid<int16_t>(5, 5, resolution, true, shift, int16_t{0});
        grid(0, 0) = -32768;
        grid(2, 2) = 0;
        grid(4, 4) = 32767;

        rastkit::RasterCollection rc;
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = resolution;

        rastkit::Layer layer;
        layer.grid = std::move(grid);
        layer.width = 5;
        layer.height = 5;
        layer.samplesPerPixel = 1;
        layer.planarConfig = 1;
        layer.datum = datum;
        layer.shift = shift;
        layer.resolution = resolution;
        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_int16.tif";
        rastkit::WriteRasterCollection(rc, testFile);

        auto rc2 = rastkit::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);
        CHECK(rc2.layers[0].holdsType<int16_t>());

        auto &g = rc2.layers[0].gridAs<int16_t>();
        CHECK(g(0, 0) == -32768);
        CHECK(g(2, 2) == 0);
        CHECK(g(4, 4) == 32767);

        std::filesystem::remove(testFile);
    }

    SUBCASE("float round-trip") {
        auto grid = dp::make_grid<float>(5, 5, resolution, true, shift, 0.0f);
        grid(0, 0) = -1.5e10f;
        grid(2, 2) = 3.14159f;
        grid(4, 4) = 1.0e-10f;

        rastkit::RasterCollection rc;
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = resolution;

        rastkit::Layer layer;
        layer.grid = std::move(grid);
        layer.width = 5;
        layer.height = 5;
        layer.samplesPerPixel = 1;
        layer.planarConfig = 1;
        layer.datum = datum;
        layer.shift = shift;
        layer.resolution = resolution;
        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_float.tif";
        rastkit::WriteRasterCollection(rc, testFile);

        auto rc2 = rastkit::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);
        CHECK(rc2.layers[0].holdsType<float>());

        auto &g = rc2.layers[0].gridAs<float>();
        CHECK(g(0, 0) == doctest::Approx(-1.5e10f));
        CHECK(g(2, 2) == doctest::Approx(3.14159f));
        CHECK(g(4, 4) == doctest::Approx(1.0e-10f));

        std::filesystem::remove(testFile);
    }

    SUBCASE("double round-trip") {
        auto grid = dp::make_grid<double>(5, 5, resolution, true, shift, 0.0);
        grid(0, 0) = -1.5e100;
        grid(2, 2) = 3.141592653589793;
        grid(4, 4) = 1.0e-100;

        rastkit::RasterCollection rc;
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = resolution;

        rastkit::Layer layer;
        layer.grid = std::move(grid);
        layer.width = 5;
        layer.height = 5;
        layer.samplesPerPixel = 1;
        layer.planarConfig = 1;
        layer.datum = datum;
        layer.shift = shift;
        layer.resolution = resolution;
        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_double.tif";
        rastkit::WriteRasterCollection(rc, testFile);

        auto rc2 = rastkit::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);
        CHECK(rc2.layers[0].holdsType<double>());

        auto &g = rc2.layers[0].gridAs<double>();
        CHECK(g(0, 0) == doctest::Approx(-1.5e100));
        CHECK(g(2, 2) == doctest::Approx(3.141592653589793));
        CHECK(g(4, 4) == doctest::Approx(1.0e-100));

        std::filesystem::remove(testFile);
    }

    SUBCASE("uint32_t round-trip") {
        auto grid = dp::make_grid<uint32_t>(5, 5, resolution, true, shift, uint32_t{0});
        grid(0, 0) = 0;
        grid(2, 2) = 2147483648U;
        grid(4, 4) = 4294967295U;

        rastkit::RasterCollection rc;
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = resolution;

        rastkit::Layer layer;
        layer.grid = std::move(grid);
        layer.width = 5;
        layer.height = 5;
        layer.samplesPerPixel = 1;
        layer.planarConfig = 1;
        layer.datum = datum;
        layer.shift = shift;
        layer.resolution = resolution;
        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_uint32.tif";
        rastkit::WriteRasterCollection(rc, testFile);

        auto rc2 = rastkit::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);
        CHECK(rc2.layers[0].holdsType<uint32_t>());

        auto &g = rc2.layers[0].gridAs<uint32_t>();
        CHECK(g(0, 0) == 0);
        CHECK(g(2, 2) == 2147483648U);
        CHECK(g(4, 4) == 4294967295U);

        std::filesystem::remove(testFile);
    }

    SUBCASE("int32_t round-trip") {
        auto grid = dp::make_grid<int32_t>(5, 5, resolution, true, shift, int32_t{0});
        grid(0, 0) = -2147483648;
        grid(2, 2) = 0;
        grid(4, 4) = 2147483647;

        rastkit::RasterCollection rc;
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = resolution;

        rastkit::Layer layer;
        layer.grid = std::move(grid);
        layer.width = 5;
        layer.height = 5;
        layer.samplesPerPixel = 1;
        layer.planarConfig = 1;
        layer.datum = datum;
        layer.shift = shift;
        layer.resolution = resolution;
        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_int32.tif";
        rastkit::WriteRasterCollection(rc, testFile);

        auto rc2 = rastkit::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);
        CHECK(rc2.layers[0].holdsType<int32_t>());

        auto &g = rc2.layers[0].gridAs<int32_t>();
        CHECK(g(0, 0) == -2147483648);
        CHECK(g(2, 2) == 0);
        CHECK(g(4, 4) == 2147483647);

        std::filesystem::remove(testFile);
    }

    SUBCASE("RGBA round-trip") {
        auto grid = dp::make_grid<rastkit::RGBA>(5, 5, resolution, true, shift, rastkit::RGBA{});
        grid(0, 0) = rastkit::RGBA{255, 0, 0, 255};     // Red
        grid(2, 2) = rastkit::RGBA{0, 255, 0, 128};     // Green with 50% alpha
        grid(4, 4) = rastkit::RGBA{0, 0, 255, 0};       // Blue fully transparent
        grid(1, 3) = rastkit::RGBA{128, 128, 128, 255}; // Gray

        rastkit::RasterCollection rc;
        rc.datum = datum;
        rc.shift = shift;
        rc.resolution = resolution;

        rastkit::Layer layer;
        layer.grid = std::move(grid);
        layer.width = 5;
        layer.height = 5;
        layer.samplesPerPixel = 4; // RGBA
        layer.planarConfig = 1;
        layer.datum = datum;
        layer.shift = shift;
        layer.resolution = resolution;
        rc.layers.push_back(std::move(layer));

        std::string testFile = "test_rgba.tif";
        rastkit::WriteRasterCollection(rc, testFile);

        auto rc2 = rastkit::ReadRasterCollection(testFile);
        REQUIRE(rc2.layers.size() == 1);
        CHECK(rc2.layers[0].holdsType<rastkit::RGBA>());
        CHECK(rc2.layers[0].isColorLayer());

        auto &g = rc2.layers[0].gridAs<rastkit::RGBA>();
        CHECK(g(0, 0).r == 255);
        CHECK(g(0, 0).g == 0);
        CHECK(g(0, 0).b == 0);
        CHECK(g(0, 0).a == 255);

        CHECK(g(2, 2).r == 0);
        CHECK(g(2, 2).g == 255);
        CHECK(g(2, 2).b == 0);
        CHECK(g(2, 2).a == 128);

        CHECK(g(4, 4).r == 0);
        CHECK(g(4, 4).g == 0);
        CHECK(g(4, 4).b == 255);
        CHECK(g(4, 4).a == 0);

        CHECK(g(1, 3).r == 128);
        CHECK(g(1, 3).g == 128);
        CHECK(g(1, 3).b == 128);
        CHECK(g(1, 3).a == 255);

        std::filesystem::remove(testFile);
    }
}

TEST_CASE("Custom tag validation") {
    SUBCASE("is_valid_custom_tag returns correct values") {
        // Reserved TIFF tags (0-32767) should be invalid
        CHECK_FALSE(rastkit::is_valid_custom_tag(0));
        CHECK_FALSE(rastkit::is_valid_custom_tag(256));   // ImageWidth
        CHECK_FALSE(rastkit::is_valid_custom_tag(339));   // SampleFormat
        CHECK_FALSE(rastkit::is_valid_custom_tag(32767)); // Last reserved

        // Private range (32768+) should be valid
        CHECK(rastkit::is_valid_custom_tag(32768)); // First private
        CHECK(rastkit::is_valid_custom_tag(50000)); // rastkit reserved start
        CHECK(rastkit::is_valid_custom_tag(50100)); // Global properties base
        CHECK(rastkit::is_valid_custom_tag(50999)); // rastkit reserved end
        CHECK(rastkit::is_valid_custom_tag(65535)); // Max tag
    }

    SUBCASE("validate_custom_tag throws for reserved tags") {
        CHECK_THROWS_AS(rastkit::validate_custom_tag(0), std::runtime_error);
        CHECK_THROWS_AS(rastkit::validate_custom_tag(256), std::runtime_error);
        CHECK_THROWS_AS(rastkit::validate_custom_tag(32767), std::runtime_error);

        CHECK_NOTHROW(rastkit::validate_custom_tag(32768));
        CHECK_NOTHROW(rastkit::validate_custom_tag(50100));
    }

    SUBCASE("Layer.setCustomTag validates tag numbers") {
        rastkit::Layer layer;

        // Valid tags should work
        CHECK_NOTHROW(layer.setCustomTag(50000, 123));
        CHECK_NOTHROW(layer.setCustomTag(32768, std::vector<uint32_t>{1, 2, 3}));

        // Invalid tags should throw
        CHECK_THROWS_AS(layer.setCustomTag(256, 123), std::runtime_error);
        CHECK_THROWS_AS(layer.setCustomTag(0, 123), std::runtime_error);
    }

    SUBCASE("Layer.getCustomTag and hasCustomTag work correctly") {
        rastkit::Layer layer;
        layer.setCustomTag(50001, 42);
        layer.setCustomTag(50002, std::vector<uint32_t>{10, 20, 30});

        CHECK(layer.hasCustomTag(50001));
        CHECK(layer.hasCustomTag(50002));
        CHECK_FALSE(layer.hasCustomTag(50003));

        auto val1 = layer.getCustomTag(50001);
        REQUIRE(val1.size() == 1);
        CHECK(val1[0] == 42);

        auto val2 = layer.getCustomTag(50002);
        REQUIRE(val2.size() == 3);
        CHECK(val2[0] == 10);
        CHECK(val2[1] == 20);
        CHECK(val2[2] == 30);

        auto val3 = layer.getCustomTag(50003);
        CHECK(val3.empty());
    }

    SUBCASE("Global properties still work (in rastkit reserved range)") {
        rastkit::Layer layer;

        // Global properties should work without throwing
        CHECK_NOTHROW(layer.setGlobalProperty("test_key", "test_value"));

        auto props = layer.getGlobalProperties();
        CHECK(props.count("test_key") == 1);
        CHECK(props["test_key"] == "test_value");
    }
}

TEST_CASE("Layer structure tests") {
    SUBCASE("Default Layer construction") {
        rastkit::Layer layer;
        CHECK(layer.ifdOffset == 0);
        CHECK(layer.width == 0);
        CHECK(layer.height == 0);
        CHECK(layer.samplesPerPixel == 0);
        CHECK(layer.planarConfig == 0);
        CHECK(layer.stripOffsets.empty());
        CHECK(layer.stripByteCounts.empty());
    }

    SUBCASE("Layer with dimensions") {
        rastkit::Layer layer;
        layer.width = 100;
        layer.height = 50;
        layer.samplesPerPixel = 1;
        layer.planarConfig = 1;

        CHECK(layer.width == 100);
        CHECK(layer.height == 50);
        CHECK(layer.samplesPerPixel == 1);
        CHECK(layer.planarConfig == 1);
    }
}

TEST_CASE("RasterCollection structure tests") {
    SUBCASE("Default RasterCollection construction") {
        rastkit::RasterCollection rc;
        CHECK(rc.layers.empty());
        // CRS is always WGS84 by default
        CHECK(rc.datum.latitude == 0.0);
        CHECK(rc.datum.longitude == 0.0);
        CHECK(rc.datum.altitude == 0.0);
        CHECK(rc.shift.rotation.to_euler().roll == 0.0);
        CHECK(rc.shift.rotation.to_euler().pitch == 0.0);
        CHECK(rc.shift.rotation.to_euler().yaw == 0.0);
    }

    SUBCASE("RasterCollection with custom values") {
        rastkit::RasterCollection rc;
        // CRS is always WGS84
        rc.datum = dp::Geo{48.0, 11.0, 500.0};
        double yaw_rad = 45.0 * M_PI / 180.0; // Convert degrees to radians
        auto rotation = dp::Quaternion::from_euler(0, 0, yaw_rad);
        rc.shift = dp::Pose{dp::Point{0, 0, 0}, rotation};
        rc.resolution = 2.0;

        // CRS is always WGS84
        CHECK(rc.datum.latitude == 48.0);
        CHECK(rc.datum.longitude == 11.0);
        CHECK(rc.datum.altitude == 500.0);
        CHECK(rc.shift.rotation.to_euler().yaw == doctest::Approx(yaw_rad).epsilon(0.001));
        CHECK(rc.resolution == 2.0);
    }

    SUBCASE("RasterCollection with layers") {
        rastkit::RasterCollection rc;

        rastkit::Layer layer1, layer2;
        layer1.width = 100;
        layer1.height = 50;
        layer2.width = 200;
        layer2.height = 100;

        rc.layers.push_back(layer1);
        rc.layers.push_back(layer2);

        CHECK(rc.layers.size() == 2);
        CHECK(rc.layers[0].width == 100);
        CHECK(rc.layers[0].height == 50);
        CHECK(rc.layers[1].width == 200);
        CHECK(rc.layers[1].height == 100);
    }
}
