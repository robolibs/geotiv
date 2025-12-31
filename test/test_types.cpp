#include <datapod/datapod.hpp>
namespace dp = datapod;
#include "geotiv/geotiv.hpp"
#include <cmath>
#include <doctest/doctest.h>

TEST_CASE("Layer structure tests") {
    SUBCASE("Default Layer construction") {
        geotiv::Layer layer;
        CHECK(layer.ifdOffset == 0);
        CHECK(layer.width == 0);
        CHECK(layer.height == 0);
        CHECK(layer.samplesPerPixel == 0);
        CHECK(layer.planarConfig == 0);
        CHECK(layer.stripOffsets.empty());
        CHECK(layer.stripByteCounts.empty());
    }

    SUBCASE("Layer with dimensions") {
        geotiv::Layer layer;
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
        geotiv::RasterCollection rc;
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
        geotiv::RasterCollection rc;
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
        geotiv::RasterCollection rc;

        geotiv::Layer layer1, layer2;
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
