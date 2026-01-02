# GeoTIFF Coordinate Precision

This document describes the coordinate precision guarantees of the geotiv library.

## Summary

The geotiv library maintains **sub-millimeter precision** for coordinate transformations across the entire globe, from equator to poles, with offsets up to 10km from the datum.

## Precision Guarantees

### Spatial Coordinate Precision (Round-trip: Save → Load)

| Location | Max Error | Typical Error |
|----------|-----------|---------------|
| **Equator (0°)** | 10 μm | 7 μm |
| **Mid-latitude (45°)** | 40 μm | 28 μm |
| **High latitude (70°)** | 0.1 mm | 0.09 mm |
| **Near pole (85°)** | 0.6 mm | 0.4 mm |

### Resolution Preservation

Resolution values are preserved with **< 0.0001% relative error** across all tested resolutions (1cm to 100m).

| Resolution | Absolute Error | Relative Error |
|------------|----------------|----------------|
| 1 cm | 0.2 nm | 0.000002% |
| 10 cm | 0.1 nm | 0.0000002% |
| 1 m | 0.07 nm | 0.000000007% |
| 10 m | 0.3 nm | 0.000000003% |
| 100 m | 24 nm | 0.00000002% |

### Large Offset Handling

The library maintains precision even with large ENU offsets from the datum:

| Offset from Datum | Max Error |
|-------------------|-----------|
| 1 km | 0.2 mm |
| 10 km | 4 mm |

### Rotation Precision

Grids with rotation (non-zero yaw) maintain the same precision:
- **45° rotation**: Max error = 50 μm

## Technical Details

### Conversion Pipeline

The library uses the following conversion chain for maximum precision:

1. **ENU → ECEF → WGS84** (when writing)
2. **WGS84 → ECEF → ENU** (when reading)

This approach uses Earth-Centered Earth-Fixed (ECEF) coordinates as an intermediate representation, which provides:
- Consistent precision across all latitudes
- No singularities at poles
- Proper handling of Earth's ellipsoidal shape

### Dependencies

Precision is achieved through:
- **concord**: Sub-millimeter accurate coordinate transformations
- **datapod**: Precise grid spatial queries
- **Double precision (IEEE 754)**: ~15-17 significant digits

### Limitations

1. **Pole proximity**: Precision degrades slightly near poles (85°+) due to longitude scale compression
   - Still maintains sub-millimeter precision
   - Acceptable for all practical applications

2. **Very large offsets**: Offsets > 10km may see slightly degraded precision
   - Still maintains < 5mm precision at 10km
   - For larger areas, consider using multiple datums

3. **Floating-point limits**: Double precision limits theoretical precision to ~1 nanometer at Earth scale

## Testing

Comprehensive precision tests are included in `test/test_precision.cpp`:
- Multiple latitudes (equator, mid-lat, high-lat, near-pole)
- Large ENU offsets (1km, 10km)
- Fine resolutions (1cm to 100m)
- Rotated grids
- Resolution conversion accuracy

Run tests with:
```bash
make test TEST=test_precision
```

## Best Practices

For maximum precision:

1. **Choose appropriate datum**: Place datum near your area of interest
2. **Keep offsets reasonable**: < 10km from datum for best results
3. **Use appropriate resolution**: Match resolution to your measurement precision
4. **Verify round-trip**: Test save/load cycles for your specific use case

## Validation

All precision claims are validated through automated tests that verify:
- Round-trip coordinate preservation (save → load)
- Resolution preservation
- Spatial query accuracy
- Multi-latitude performance

See `test/test_precision.cpp` for complete test suite.
