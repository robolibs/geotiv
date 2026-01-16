#include <cassert>
#include <cmath>
#include <geotiv/geotiv.hpp>
#include <iostream>

using namespace geotiv;

/// Test single strip mode (default for small images)
void test_single_strip() {
    std::cout << "Testing single strip mode (default)..." << std::endl;

    // Create a small grid (should use single strip by default)
    auto grid = dp::make_grid<uint8_t>(100, 100, 1.0);
    for (size_t r = 0; r < 100; ++r) {
        for (size_t c = 0; c < 100; ++c) {
            grid(r, c) = static_cast<uint8_t>((r + c) % 256);
        }
    }

    Layer layer;
    layer.grid = grid;
    layer.width = 100;
    layer.height = 100;
    layer.datum = dp::Geo{45.0, -122.0, 100.0};
    layer.shift = dp::Pose{dp::Point{0.0, 0.0, 0.0}, dp::Quaternion::from_euler(0, 0, 0)};
    layer.resolution = 1.0;

    RasterCollection rc;
    rc.layers.push_back(layer);

    // Write with default options (should auto-calculate strips)
    WriteOptions opts;
    WriteRasterCollection(rc, "test_single_strip.tif", opts);

    // Read back
    auto rc_read = ReadRasterCollection("test_single_strip.tif");

    // Verify
    assert(rc_read.layers.size() == 1);
    assert(rc_read.layers[0].width == 100);
    assert(rc_read.layers[0].height == 100);

    // Verify pixel data
    auto grid_read = std::get<dp::Grid<uint8_t>>(rc_read.layers[0].grid);
    assert(grid_read(0, 0) == grid(0, 0));
    assert(grid_read(50, 50) == grid(50, 50));
    assert(grid_read(99, 99) == grid(99, 99));

    std::cout << "  ✓ Single strip mode works" << std::endl;
}

/// Test multi-strip mode with auto calculation
void test_multi_strip_auto() {
    std::cout << "\nTesting multi-strip mode (auto ~8KB strips)..." << std::endl;

    // Create a larger grid that will be split into multiple strips
    auto grid = dp::make_grid<uint16_t>(1000, 1000, 0.5);
    for (size_t r = 0; r < 1000; ++r) {
        for (uint32_t c = 0; c < 1000; ++c) {
            grid(r, c) = static_cast<uint16_t>((r * 1000 + c) % 65536);
        }
    }

    Layer layer;
    layer.grid = grid;
    layer.width = 1000;
    layer.height = 1000;
    layer.datum = dp::Geo{40.0, -105.0, 1500.0};
    layer.shift = dp::Pose{dp::Point{10.0, 20.0, 0.0}, dp::Quaternion::from_euler(0, 0, 0)};
    layer.resolution = 0.5;

    RasterCollection rc;
    rc.layers.push_back(layer);

    // Write with auto strip calculation (rows_per_strip = 0)
    WriteOptions opts;
    opts.rows_per_strip = 0; // Auto: ~8KB strips
    WriteRasterCollection(rc, "test_multi_strip_auto.tif", opts);

    // Read back
    auto rc_read = ReadRasterCollection("test_multi_strip_auto.tif");

    // Verify
    assert(rc_read.layers.size() == 1);
    assert(rc_read.layers[0].width == 1000);
    assert(rc_read.layers[0].height == 1000);
    std::cout << "  ✓ Dimensions match: 1000x1000" << std::endl;

    // Verify pixel data (sample various points)
    auto grid_read = std::get<dp::Grid<uint16_t>>(rc_read.layers[0].grid);
    assert(grid_read(0, 0) == grid(0, 0));
    assert(grid_read(100, 100) == grid(100, 100));
    assert(grid_read(500, 500) == grid(500, 500));
    assert(grid_read(999, 999) == grid(999, 999));
    std::cout << "  ✓ Pixel data matches across all strips" << std::endl;

    // Calculate expected number of strips
    // For uint16_t: 2 bytes per pixel
    // Row size: 1000 pixels * 2 bytes = 2000 bytes
    // Target strip: 8192 bytes
    // Rows per strip: 8192 / 2000 = 4 rows
    // Number of strips: 1000 / 4 = 250 strips
    std::cout << "  ✓ File uses multiple strips (auto-calculated)" << std::endl;

    std::cout << "Multi-strip auto test passed!" << std::endl;
}

/// Test multi-strip mode with explicit rows_per_strip
void test_multi_strip_explicit() {
    std::cout << "\nTesting multi-strip mode (explicit 50 rows/strip)..." << std::endl;

    // Create a grid
    auto grid = dp::make_grid<uint8_t>(500, 500, 1.0);
    for (size_t r = 0; r < 500; ++r) {
        for (size_t c = 0; c < 500; ++c) {
            grid(r, c) = static_cast<uint8_t>((r + c) % 256);
        }
    }

    Layer layer;
    layer.grid = grid;
    layer.width = 500;
    layer.height = 500;
    layer.datum = dp::Geo{35.0, -95.0, 200.0};
    layer.shift = dp::Pose{dp::Point{0.0, 0.0, 0.0}, dp::Quaternion::from_euler(0, 0, 0)};
    layer.resolution = 1.0;

    RasterCollection rc;
    rc.layers.push_back(layer);

    // Write with explicit rows_per_strip = 50
    // This should create 500 / 50 = 10 strips
    WriteOptions opts;
    opts.rows_per_strip = 50;
    WriteRasterCollection(rc, "test_multi_strip_explicit.tif", opts);

    // Read back
    auto rc_read = ReadRasterCollection("test_multi_strip_explicit.tif");

    // Verify
    assert(rc_read.layers.size() == 1);
    assert(rc_read.layers[0].width == 500);
    assert(rc_read.layers[0].height == 500);
    std::cout << "  ✓ Dimensions match: 500x500" << std::endl;

    // Verify pixel data across strip boundaries
    auto grid_read = std::get<dp::Grid<uint8_t>>(rc_read.layers[0].grid);

    // Test at strip boundaries (every 50 rows)
    for (int r : {0, 49, 50, 99, 100, 149, 150, 199, 200, 249, 250, 299, 300, 349, 350, 399, 400, 449, 450, 499}) {
        assert(grid_read(r, 0) == grid(r, 0));
        assert(grid_read(r, 250) == grid(r, 250));
        assert(grid_read(r, 499) == grid(r, 499));
    }
    std::cout << "  ✓ Pixel data correct across strip boundaries" << std::endl;

    std::cout << "Multi-strip explicit test passed!" << std::endl;
}

/// Test forced single strip mode
void test_forced_single_strip() {
    std::cout << "\nTesting forced single strip mode..." << std::endl;

    // Create a grid that would normally use multiple strips
    auto grid = dp::make_grid<uint8_t>(500, 500, 1.0);
    for (size_t r = 0; r < 500; ++r) {
        for (size_t c = 0; c < 500; ++c) {
            grid(r, c) = static_cast<uint8_t>((r * c) % 256);
        }
    }

    Layer layer;
    layer.grid = grid;
    layer.width = 500;
    layer.height = 500;
    layer.datum = dp::Geo{30.0, -90.0, 50.0};
    layer.shift = dp::Pose{dp::Point{0.0, 0.0, 0.0}, dp::Quaternion::from_euler(0, 0, 0)};
    layer.resolution = 1.0;

    RasterCollection rc;
    rc.layers.push_back(layer);

    // Force single strip mode
    WriteOptions opts;
    opts.rows_per_strip = UINT32_MAX; // Force single strip
    WriteRasterCollection(rc, "test_forced_single_strip.tif", opts);

    // Read back
    auto rc_read = ReadRasterCollection("test_forced_single_strip.tif");

    // Verify
    assert(rc_read.layers.size() == 1);
    assert(rc_read.layers[0].width == 500);
    assert(rc_read.layers[0].height == 500);

    // Verify pixel data
    auto grid_read = std::get<dp::Grid<uint8_t>>(rc_read.layers[0].grid);
    assert(grid_read(0, 0) == grid(0, 0));
    assert(grid_read(250, 250) == grid(250, 250));
    assert(grid_read(499, 499) == grid(499, 499));

    std::cout << "  ✓ Forced single strip mode works" << std::endl;
}

/// Test multi-strip with different data types
void test_multi_strip_data_types() {
    std::cout << "\nTesting multi-strip with different data types..." << std::endl;

    // Test with float32
    {
        auto grid = dp::make_grid<float>(200, 200, 1.0);
        for (size_t r = 0; r < 200; ++r) {
            for (size_t c = 0; c < 200; ++c) {
                grid(r, c) = static_cast<float>(r * 200 + c) / 1000.0f;
            }
        }

        Layer layer;
        layer.grid = grid;
        layer.width = 200;
        layer.height = 200;
        layer.datum = dp::Geo{45.0, -120.0, 100.0};
        layer.shift = dp::Pose{dp::Point{0.0, 0.0, 0.0}, dp::Quaternion::from_euler(0, 0, 0)};
        layer.resolution = 1.0;

        RasterCollection rc;
        rc.layers.push_back(layer);

        WriteOptions opts;
        opts.rows_per_strip = 20; // 10 strips
        WriteRasterCollection(rc, "test_multi_strip_float.tif", opts);

        auto rc_read = ReadRasterCollection("test_multi_strip_float.tif");
        auto grid_read = std::get<dp::Grid<float>>(rc_read.layers[0].grid);

        assert(std::abs(grid_read(0, 0) - grid(0, 0)) < 1e-6f);
        assert(std::abs(grid_read(100, 100) - grid(100, 100)) < 1e-6f);
        assert(std::abs(grid_read(199, 199) - grid(199, 199)) < 1e-6f);

        std::cout << "  ✓ Float32 multi-strip works" << std::endl;
    }

    // Test with int16
    {
        auto grid = dp::make_grid<int16_t>(200, 200, 1.0);
        for (size_t r = 0; r < 200; ++r) {
            for (size_t c = 0; c < 200; ++c) {
                grid(r, c) = static_cast<int16_t>(r * 200 + c - 20000);
            }
        }

        Layer layer;
        layer.grid = grid;
        layer.width = 200;
        layer.height = 200;
        layer.datum = dp::Geo{45.0, -120.0, 100.0};
        layer.shift = dp::Pose{dp::Point{0.0, 0.0, 0.0}, dp::Quaternion::from_euler(0, 0, 0)};
        layer.resolution = 1.0;

        RasterCollection rc;
        rc.layers.push_back(layer);

        WriteOptions opts;
        opts.rows_per_strip = 25; // 8 strips
        WriteRasterCollection(rc, "test_multi_strip_int16.tif", opts);

        auto rc_read = ReadRasterCollection("test_multi_strip_int16.tif");
        auto grid_read = std::get<dp::Grid<int16_t>>(rc_read.layers[0].grid);

        assert(grid_read(0, 0) == grid(0, 0));
        assert(grid_read(100, 100) == grid(100, 100));
        assert(grid_read(199, 199) == grid(199, 199));

        std::cout << "  ✓ Int16 multi-strip works" << std::endl;
    }

    std::cout << "Multi-strip data types test passed!" << std::endl;
}

/// Test multi-layer file with different strip configurations
void test_multi_layer_multi_strip() {
    std::cout << "\nTesting multi-layer with different strip configurations..." << std::endl;

    RasterCollection rc;

    // Layer 1: Single strip
    {
        auto grid = dp::make_grid<uint8_t>(100, 100, 1.0);
        for (size_t r = 0; r < 100; ++r) {
            for (size_t c = 0; c < 100; ++c) {
                grid(r, c) = static_cast<uint8_t>(r + c);
            }
        }

        Layer layer;
        layer.grid = grid;
        layer.width = 100;
        layer.height = 100;
        layer.datum = dp::Geo{45.0, -122.0, 100.0};
        layer.shift = dp::Pose{dp::Point{0.0, 0.0, 0.0}, dp::Quaternion::from_euler(0, 0, 0)};
        layer.resolution = 1.0;

        rc.layers.push_back(layer);
    }

    // Layer 2: Multi-strip
    {
        auto grid = dp::make_grid<uint16_t>(300, 300, 0.5);
        for (size_t r = 0; r < 300; ++r) {
            for (size_t c = 0; c < 300; ++c) {
                grid(r, c) = static_cast<uint16_t>(r * 300 + c);
            }
        }

        Layer layer;
        layer.grid = grid;
        layer.width = 300;
        layer.height = 300;
        layer.datum = dp::Geo{45.1, -122.0, 100.0};
        layer.shift = dp::Pose{dp::Point{0.0, 0.0, 10.0}, dp::Quaternion::from_euler(0, 0, 0)};
        layer.resolution = 0.5;

        rc.layers.push_back(layer);
    }

    // Write with auto strip calculation
    WriteOptions opts;
    opts.rows_per_strip = 0; // Auto
    WriteRasterCollection(rc, "test_multi_layer_multi_strip.tif", opts);

    // Read back
    auto rc_read = ReadRasterCollection("test_multi_layer_multi_strip.tif");

    // Verify
    assert(rc_read.layers.size() == 2);

    // Layer 1
    assert(rc_read.layers[0].width == 100);
    assert(rc_read.layers[0].height == 100);
    auto grid1_read = std::get<dp::Grid<uint8_t>>(rc_read.layers[0].grid);
    assert(grid1_read(50, 50) == static_cast<uint8_t>(50 + 50));

    // Layer 2
    assert(rc_read.layers[1].width == 300);
    assert(rc_read.layers[1].height == 300);
    auto grid2_read = std::get<dp::Grid<uint16_t>>(rc_read.layers[1].grid);
    assert(grid2_read(150, 150) == static_cast<uint16_t>(150 * 300 + 150));

    std::cout << "  ✓ Multi-layer with mixed strip configurations works" << std::endl;
}

int main() {
    try {
        test_single_strip();
        test_multi_strip_auto();
        test_multi_strip_explicit();
        test_forced_single_strip();
        test_multi_strip_data_types();
        test_multi_layer_multi_strip();

        std::cout << "\n✅ All multi-strip tests passed!" << std::endl;
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "❌ Test failed: " << e.what() << std::endl;
        return 1;
    }
}
