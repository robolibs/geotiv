# TODO for Next Session

## Priority 1: Complete Metadata Tags (geotiv-2vu.1, 2vu.2, 2vu.3)

These were implemented and tested but lost during git operations. Need to re-implement:

### 1. Software Tag (305)
- [x] Add to WriteOptions ✅
- [ ] Add softwareOffsets vector
- [ ] Calculate offset in variable-length data section
- [ ] Add IFD entry: `{305, 2, softwareLength, softwareOffsets[i]}`
- [ ] Write software string data
- [ ] Test: verify tag 305 exists in output

### 2. DateTime Tag (306)
- [x] Add to WriteOptions ✅
- [x] Add get_current_datetime() helper ✅
- [ ] Add datetimeOffsets vector
- [ ] Calculate offset in variable-length data section
- [ ] Add IFD entry: `{306, 2, datetimeLength, datetimeOffsets[i]}`
- [ ] Write datetime string data
- [ ] Test: verify tag 306 exists and format is correct

### 3. XResolution/YResolution Tags (282/283/296)
- [x] Add to WriteOptions ✅
- [ ] Add xresolutionOffsets, yresolutionOffsets vectors
- [ ] Calculate offsets (8 bytes each for RATIONAL type)
- [ ] Add IFD entries:
  - `{282, 5, 1, xresolutionOffsets[i]}` (RATIONAL)
  - `{283, 5, 1, yresolutionOffsets[i]}` (RATIONAL)
  - `{296, 3, 1, options.resolution_unit}` (SHORT)
- [ ] Write RATIONAL data (numerator, denominator pairs)
- [ ] Test: verify tags exist and values are correct

**Reference**: See commit history or test_writer.cpp for test cases that were written.

## Priority 2: Multi-Strip Organization (geotiv-2vu.4)

### Writer Changes:

1. **Calculate optimal strip size**:
   ```cpp
   uint32_t calculate_rows_per_strip(uint32_t width, uint32_t bytes_per_pixel) {
       const uint32_t target_strip_bytes = 8192; // 8KB
       uint32_t bytes_per_row = width * bytes_per_pixel;
       uint32_t rows = target_strip_bytes / bytes_per_row;
       return std::max(1u, rows); // At least 1 row
   }
   ```

2. **Split data into strips**:
   - Change from `std::vector<std::vector<uint8_t>> strips(N)` (one per layer)
   - To: `std::vector<std::vector<std::vector<uint8_t>>> layerStrips(N)` (multiple per layer)
   - Split grid data row-by-row into strips

3. **Update IFD entries**:
   - StripOffsets: array of offsets (one per strip)
   - StripByteCounts: array of byte counts (one per strip)
   - RowsPerStrip: calculated value (not image height)

### Parser Changes:

Current parser already handles multiple strips! Just verify it works:
```cpp
for (size_t i = 0; i < L.stripOffsets.size(); ++i) {
    f.seekg(L.stripOffsets[i], std::ios::beg);
    f.read(...);
}
```

### Testing:
- Create large image (e.g., 1000x1000)
- Verify multiple strips created
- Verify round-trip works
- Check memory usage is lower

## Priority 3: ROI (Region of Interest) Loading

Add new function to read only a region:

```cpp
struct ROI {
    uint32_t row_start;
    uint32_t row_end;
    uint32_t col_start;
    uint32_t col_end;
};

RasterCollection ReadRasterRegion(const fs::path &file, const ROI &roi);
```

Implementation:
1. Read TIFF header and IFD
2. Calculate which strips overlap with ROI
3. Read only those strips
4. Extract ROI data from strips
5. Create smaller grid with correct spatial info

**Example Use Case**:
```cpp
// Load 500x500 region around point (5000, 3000) from 10000x10000 image
ROI roi{2750, 3250, 4750, 5250};
auto rc = ReadRasterRegion("huge_image.tif", roi);
// Only reads ~40 strips instead of entire file!
```

## Files to Modify:

- `include/geotiv/writter.hpp` - metadata tags, multi-strip writing
- `include/geotiv/parser.hpp` - ROI loading function
- `test/test_writer.cpp` - metadata tag tests
- `test/test_multi_strip.cpp` - NEW: multi-strip tests
- `test/test_roi.cpp` - NEW: ROI loading tests

## Estimated Time:

- Metadata tags: 30-45 minutes (straightforward, already done once)
- Multi-strip: 1-2 hours (refactoring data structures)
- ROI loading: 1 hour (new feature)
- Testing: 30 minutes

**Total: 3-4 hours**

## Notes:

- All metadata tag tests were written and working before
- Parser already supports multi-strip reading
- ROI loading is the most complex new feature
- Consider adding progress callback for large file operations
