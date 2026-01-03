# TODO for Next Session

## Session Summary (2026-01-03)

### ✅ Completed This Session:

1. **BigTIFF Support (geotiv-2vu.5)** - COMPLETE
   - ✅ Parser: Detect magic 43, read 64-bit offsets
   - ✅ Writer: Auto-detect >4GB or force with `WriteOptions.force_bigtiff`
   - ✅ BigTIFF header: 16 bytes with 64-bit IFD offset
   - ✅ BigTIFF IFD: 64-bit entry count, 20-byte entries
   - ✅ Comprehensive tests: header, round-trip, multi-layer
   - ✅ Verified with gdalinfo
   - ✅ All 12 tests passing

2. **Enhanced Georeferencing Epic (geotiv-qya)** - CLOSED
   - ✅ All 5 tasks complete (100%)
   - ✅ GeoDoubleParamsTag (34736)
   - ✅ GeoAsciiParamsTag (34737)
   - ✅ Vertical CRS GeoKeys
   - ✅ Citation GeoKeys
   - ✅ Coordinate precision handling

### 📊 Project Status:

**Main Epic (geotiv-blp):** 5/6 child epics complete (83%)
- ✅ Critical Bug Fixes (geotiv-b2z) - CLOSED
- ✅ Extended Bit Depth (geotiv-101) - CLOSED
- ✅ Compression Support (geotiv-nxd) - CLOSED
- ✅ Color and Multi-Band (geotiv-2n1) - CLOSED
- ✅ Enhanced Georeferencing (geotiv-qya) - CLOSED
- ⏳ TIFF Metadata (geotiv-2vu) - 6/7 tasks (86%)

**Metadata Epic (geotiv-2vu):** 6/7 tasks complete (86%)
- ✅ Software tag (305)
- ✅ DateTime tag (306)
- ✅ XResolution/YResolution tags (282/283/296)
- ✅ NoData value (42113)
- ✅ Error messages infrastructure
- ✅ BigTIFF support
- ⏳ Multi-strip organization (DEFERRED - complex, 2-3 hours)

## Priority 1: Multi-Strip Organization (geotiv-2vu.4)

**Status:** Open, deferred from previous sessions
**Complexity:** High (2-3 hours)
**Priority:** P2 (not blocking other work)

### Why Multi-Strip?

Current implementation uses single strip (entire image in memory):
- Large images cause memory pressure
- No streaming read possible
- No region-of-interest (ROI) loading

Multi-strip benefits:
- Process image in chunks (~8KB strips)
- Lower memory footprint
- Enable ROI loading (load 500x500 from 10000x10000 without loading entire file)

### Implementation Plan:

#### 1. Writer Changes (`include/geotiv/writter.hpp`):

**a) Calculate optimal strip size:**
```cpp
uint32_t calculate_rows_per_strip(uint32_t width, uint32_t bytes_per_pixel, uint32_t user_rows) {
    if (user_rows == UINT32_MAX) return height; // Single strip
    if (user_rows > 0) return user_rows;        // User specified
    
    // Auto: target ~8KB strips
    const uint32_t target_strip_bytes = 8192;
    uint32_t bytes_per_row = width * bytes_per_pixel;
    uint32_t rows = target_strip_bytes / bytes_per_row;
    return std::max(1u, rows);
}
```

**b) Refactor data structures:**
- Change: `std::vector<std::vector<uint8_t>> strips(N)` (one per layer)
- To: `std::vector<std::vector<std::vector<uint8_t>>> layerStrips(N)` (multiple per layer)
- Update: `stripOffsets`, `stripCounts` to be vectors of vectors

**c) Split grid data into strips:**
```cpp
for (uint32_t strip_idx = 0; strip_idx < num_strips; ++strip_idx) {
    uint32_t start_row = strip_idx * rows_per_strip;
    uint32_t end_row = std::min(start_row + rows_per_strip, height);
    // Copy rows [start_row, end_row) into strip
}
```

**d) Update IFD entries:**
- StripOffsets: type 4 (LONG), count = num_strips, offset to array
- StripByteCounts: type 4 (LONG), count = num_strips, offset to array
- RowsPerStrip: type 4 (LONG), count = 1, value = rows_per_strip

#### 2. Parser Changes (`include/geotiv/parser.hpp`):

Parser already handles multiple strips! Just verify:
```cpp
for (size_t i = 0; i < L.stripOffsets.size(); ++i) {
    f.seekg(L.stripOffsets[i], std::ios::beg);
    f.read(reinterpret_cast<char*>(&rawData[offset]), L.stripByteCounts[i]);
    offset += L.stripByteCounts[i];
}
```

#### 3. Add ROI Loading Function:

```cpp
struct ROI {
    uint32_t row_start, row_end;
    uint32_t col_start, col_end;
};

RasterCollection ReadRasterRegion(const fs::path &file, const ROI &roi, size_t layer_idx = 0);
```

Implementation:
1. Read TIFF header and IFD for specified layer
2. Calculate which strips overlap with ROI rows
3. Read only those strips
4. Extract ROI columns from strip data
5. Create smaller grid with correct spatial info (adjust shift)

#### 4. Testing:

Create `test/test_multi_strip.cpp`:
- Test with large image (1000x1000)
- Verify multiple strips created
- Verify round-trip works
- Test ROI loading (load 100x100 from 1000x1000)
- Verify memory efficiency

### Files to Modify:
- `include/geotiv/writter.hpp` - multi-strip writing
- `include/geotiv/parser.hpp` - ROI loading function
- `test/test_multi_strip.cpp` - NEW: comprehensive tests

### Estimated Time: 2-3 hours

## Priority 2: Close Remaining Epics

Once multi-strip is complete:
1. Close geotiv-2vu epic (7/7 tasks)
2. Close geotiv-blp main epic (6/6 child epics)
3. Project reaches 100% completion! 🎉

## Priority 3: Documentation (Optional)

If time permits:
- Update README.md with BigTIFF support
- Document ROI loading feature
- Add examples for multi-strip usage
- Performance benchmarks (memory usage comparison)

## Notes:

- All previous work (metadata tags, NoData, BigTIFF, georeferencing) is complete and tested
- Multi-strip is the last remaining feature
- This is a significant refactoring but well-defined
- Parser already supports reading multi-strip files
- ROI loading is the most valuable new feature

## Quick Reference:

```bash
# Build and test
make build
make test

# Run specific test
make test TEST=test_multi_strip

# Check beads status
bd show geotiv-2vu.4
bd ready

# When complete
bd close geotiv-2vu.4 -r "Implemented multi-strip organization..."
bd close geotiv-2vu -r "All 7 tasks complete..."
bd close geotiv-blp -r "All 6 child epics complete..."
git add . && git commit -m "feat: add multi-strip organization..."
git push && bd sync
```
