#include <cassert>
#include <cmath>
#include <fstream>
#include <geotiv/geotiv.hpp>
#include <iostream>

using namespace geotiv;

/// Helper to read bytes from a file
std::vector<uint8_t> read_file_bytes(const std::string &path) {
    std::ifstream file(path, std::ios::binary | std::ios::ate);
    if (!file) {
        throw std::runtime_error("Cannot open file: " + path);
    }
    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);
    std::vector<uint8_t> buffer(size);
    if (!file.read(reinterpret_cast<char *>(buffer.data()), size)) {
        throw std::runtime_error("Cannot read file: " + path);
    }
    return buffer;
}

/// Helper to read little-endian uint16
uint16_t read_le16(const std::vector<uint8_t> &buf, size_t pos) {
    return static_cast<uint16_t>(buf[pos]) | (static_cast<uint16_t>(buf[pos + 1]) << 8);
}

/// Helper to read little-endian uint32
uint32_t read_le32(const std::vector<uint8_t> &buf, size_t pos) {
    return static_cast<uint32_t>(buf[pos]) | (static_cast<uint32_t>(buf[pos + 1]) << 8) |
           (static_cast<uint32_t>(buf[pos + 2]) << 16) | (static_cast<uint32_t>(buf[pos + 3]) << 24);
}

/// Helper to read little-endian uint64
uint64_t read_le64(const std::vector<uint8_t> &buf, size_t pos) {
    uint64_t result = 0;
    for (int i = 0; i < 8; ++i) {
        result |= (static_cast<uint64_t>(buf[pos + i]) << (8 * i));
    }
    return result;
}

/// Test BigTIFF header format
void test_bigtiff_header() {
    std::cout << "Testing BigTIFF header format..." << std::endl;

    // Create a simple grid
    auto grid = dp::make_grid<uint8_t>(100, 100, 1.0);
    for (size_t r = 0; r < 100; ++r) {
        for (size_t c = 0; c < 100; ++c) {
            grid(r, c) = static_cast<uint8_t>((r + c) % 256);
        }
    }

    // Create layer
    Layer layer;
    layer.grid = grid;
    layer.width = 100;
    layer.height = 100;
    layer.datum = dp::Geo{45.0, -122.0, 100.0};
    layer.shift = dp::Pose{dp::Point{0.0, 0.0, 0.0}, dp::Quaternion::from_euler(0, 0, 0)};
    layer.resolution = 1.0;

    RasterCollection rc;
    rc.layers.push_back(layer);

    // Write with force_bigtiff = true
    WriteOptions opts;
    opts.force_bigtiff = true;
    WriteRasterCollection(rc, "test_bigtiff_forced.tif", opts);

    // Read the file and verify BigTIFF header
    auto bytes = read_file_bytes("test_bigtiff_forced.tif");

    // Check byte order (II = little-endian)
    assert(bytes[0] == 'I' && bytes[1] == 'I');

    // Check magic number (43 for BigTIFF)
    uint16_t magic = read_le16(bytes, 2);
    assert(magic == 43);
    std::cout << "  ✓ Magic number is 43 (BigTIFF)" << std::endl;

    // Check offset size (should be 8)
    uint16_t offset_size = read_le16(bytes, 4);
    assert(offset_size == 8);
    std::cout << "  ✓ Offset size is 8 bytes" << std::endl;

    // Check reserved field (should be 0)
    uint16_t reserved = read_le16(bytes, 6);
    assert(reserved == 0);
    std::cout << "  ✓ Reserved field is 0" << std::endl;

    // Check first IFD offset (64-bit)
    uint64_t first_ifd = read_le64(bytes, 8);
    assert(first_ifd > 16); // Should be after header + pixel data
    std::cout << "  ✓ First IFD offset: " << first_ifd << " (64-bit)" << std::endl;

    std::cout << "BigTIFF header test passed!" << std::endl;
}

/// Test that classic TIFF is used by default for small files
void test_classic_tiff_default() {
    std::cout << "\nTesting classic TIFF for small files..." << std::endl;

    // Create a small grid
    auto grid = dp::make_grid<uint8_t>(50, 50, 1.0);
    for (size_t r = 0; r < 50; ++r) {
        for (size_t c = 0; c < 50; ++c) {
            grid(r, c) = static_cast<uint8_t>((r * c) % 256);
        }
    }

    Layer layer;
    layer.grid = grid;
    layer.width = 50;
    layer.height = 50;
    layer.datum = dp::Geo{45.0, -122.0, 100.0};
    layer.shift = dp::Pose{dp::Point{0.0, 0.0, 0.0}, dp::Quaternion::from_euler(0, 0, 0)};
    layer.resolution = 1.0;

    RasterCollection rc;
    rc.layers.push_back(layer);

    // Write without forcing BigTIFF
    WriteOptions opts;
    opts.force_bigtiff = false;
    WriteRasterCollection(rc, "test_classic_tiff.tif", opts);

    // Read the file and verify classic TIFF header
    auto bytes = read_file_bytes("test_classic_tiff.tif");

    // Check byte order
    assert(bytes[0] == 'I' && bytes[1] == 'I');

    // Check magic number (42 for classic TIFF)
    uint16_t magic = read_le16(bytes, 2);
    assert(magic == 42);
    std::cout << "  ✓ Magic number is 42 (classic TIFF)" << std::endl;

    // Check first IFD offset (32-bit)
    uint32_t first_ifd = read_le32(bytes, 4);
    assert(first_ifd > 8); // Should be after header + pixel data
    std::cout << "  ✓ First IFD offset: " << first_ifd << " (32-bit)" << std::endl;

    std::cout << "Classic TIFF test passed!" << std::endl;
}

/// Test BigTIFF round-trip (write and read back)
void test_bigtiff_roundtrip() {
    std::cout << "\nTesting BigTIFF round-trip..." << std::endl;

    // Create a grid with known pattern
    auto grid = dp::make_grid<uint16_t>(200, 200, 0.5);
    for (size_t r = 0; r < 200; ++r) {
        for (size_t c = 0; c < 200; ++c) {
            grid(r, c) = static_cast<uint16_t>((r * 200 + c) % 65536);
        }
    }

    Layer layer;
    layer.grid = grid;
    layer.width = 200;
    layer.height = 200;
    layer.datum = dp::Geo{40.0, -105.0, 1500.0};
    layer.shift = dp::Pose{dp::Point{10.0, 20.0, 0.0}, dp::Quaternion::from_euler(0, 0, 0)};
    layer.resolution = 0.5;
    layer.noDataValue = 65535.0;

    RasterCollection rc;
    rc.layers.push_back(layer);

    // Write as BigTIFF
    WriteOptions opts;
    opts.force_bigtiff = true;
    WriteRasterCollection(rc, "test_bigtiff_roundtrip.tif", opts);

    // Read back
    auto rc_read = ReadRasterCollection("test_bigtiff_roundtrip.tif");

    // Verify we got the data back
    assert(rc_read.layers.size() == 1);
    std::cout << "  ✓ Read back 1 layer" << std::endl;

    auto &layer_read = rc_read.layers[0];
    assert(layer_read.width == 200);
    assert(layer_read.height == 200);
    std::cout << "  ✓ Dimensions match: 200x200" << std::endl;

    // Verify datum
    assert(std::abs(layer_read.datum.latitude - 40.0) < 1e-6);
    assert(std::abs(layer_read.datum.longitude - (-105.0)) < 1e-6);
    assert(std::abs(layer_read.datum.altitude - 1500.0) < 1e-6);
    std::cout << "  ✓ Datum matches" << std::endl;

    // Verify NoData value
    assert(layer_read.noDataValue.has_value());
    assert(std::abs(layer_read.noDataValue.value() - 65535.0) < 1e-6);
    std::cout << "  ✓ NoData value matches" << std::endl;

    // Verify pixel data (sample a few points)
    auto grid_read = std::get<dp::Grid<uint16_t>>(layer_read.grid);
    assert(grid_read(0, 0) == 0);
    assert(grid_read(10, 10) == (10 * 200 + 10) % 65536);
    assert(grid_read(199, 199) == (199 * 200 + 199) % 65536);
    std::cout << "  ✓ Pixel data matches" << std::endl;

    std::cout << "BigTIFF round-trip test passed!" << std::endl;
}

/// Test multi-layer BigTIFF
void test_bigtiff_multi_layer() {
    std::cout << "\nTesting multi-layer BigTIFF..." << std::endl;

    RasterCollection rc;

    // Create 3 layers with different data
    for (int layer_idx = 0; layer_idx < 3; ++layer_idx) {
        auto grid = dp::make_grid<uint8_t>(100, 100, 1.0);
        for (size_t r = 0; r < 100; ++r) {
            for (size_t c = 0; c < 100; ++c) {
                grid(r, c) = static_cast<uint8_t>((r + c + layer_idx * 50) % 256);
            }
        }

        Layer layer;
        layer.grid = grid;
        layer.width = 100;
        layer.height = 100;
        layer.datum = dp::Geo{45.0 + layer_idx, -122.0, 100.0};
        layer.shift = dp::Pose{dp::Point{0.0, 0.0, layer_idx * 10.0}, dp::Quaternion::from_euler(0, 0, 0)};
        layer.resolution = 1.0;

        rc.layers.push_back(layer);
    }

    // Write as BigTIFF
    WriteOptions opts;
    opts.force_bigtiff = true;
    WriteRasterCollection(rc, "test_bigtiff_multi.tif", opts);

    // Read back
    auto rc_read = ReadRasterCollection("test_bigtiff_multi.tif");

    // Verify
    assert(rc_read.layers.size() == 3);
    std::cout << "  ✓ Read back 3 layers" << std::endl;

    for (int i = 0; i < 3; ++i) {
        assert(rc_read.layers[i].width == 100);
        assert(rc_read.layers[i].height == 100);
        assert(std::abs(rc_read.layers[i].datum.latitude - (45.0 + i)) < 1e-6);
    }
    std::cout << "  ✓ All layers have correct metadata" << std::endl;

    std::cout << "Multi-layer BigTIFF test passed!" << std::endl;
}

int main() {
    try {
        test_bigtiff_header();
        test_classic_tiff_default();
        test_bigtiff_roundtrip();
        test_bigtiff_multi_layer();

        std::cout << "\n✅ All BigTIFF tests passed!" << std::endl;
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "❌ Test failed: " << e.what() << std::endl;
        return 1;
    }
}
