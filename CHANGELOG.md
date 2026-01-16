# Changelog

## [0.0.6] - 2026-01-16

### <!-- 0 -->⛰️  Features

- Add convenience helpers for GridVariant

## [0.0.5] - 2026-01-16

### <!-- 0 -->⛰️  Features

- Add GeoTIFF enhancements for metadata and data types
- Improve TIFF writing and GeoTIFF support
- Feat: Bootstrap geotiv project infrastructure

### <!-- 7 -->⚙️ Miscellaneous Tasks

- Remove old TIFF test images
- Update dependencies and Makefile options

## [0.0.4] - 2026-01-16

### <!-- 7 -->⚙️ Miscellaneous Tasks

- Update dependency versions

## [0.0.3] - 2026-01-03

### <!-- 0 -->⛰️  Features

- Add multi-strip organization for memory-efficient large file handling
- Add BigTIFF support for files >4GB
- Add custom error types for better error messages
- Add vertical CRS GeoKeys support
- Add GeoDoubleParamsTag (34736) support
- Add GeoAsciiParamsTag (34737) and citation GeoKeys
- Add NoData value support (GDAL_NODATA tag 42113)
- Implement TIFF metadata tags (Software, DateTime, Resolution)
- Add comprehensive coordinate precision validation
- Support uncompressed TIFF writing
- Add RGBA color image support with pigment integration
- Add multi-type grid support with GridVariant

### <!-- 1 -->🐛 Bug Fixes

- Critical TIFF spec compliance and bug fixes

### <!-- 3 -->📚 Documentation

- Update TODO with BigTIFF completion and session summary
- Add detailed TODO for next session

### <!-- 7 -->⚙️ Miscellaneous Tasks

- Update library dependency versions
- Remove outdated files and tests

### Wip

- Add metadata tags structure and multi-strip foundation

## [0.0.2] - 2025-12-31

### <!-- 0 -->⛰️  Features

- Introduce comprehensive build system and dependency management
- Convert geotiv to a header-only library

### <!-- 2 -->🚜 Refactor

- Migrate to new concord ecosystem (datapod, optinum, concord)

### <!-- 7 -->⚙️ Miscellaneous Tasks

- Update dependencies and project version
- Clean up devbox and doctest configurations

### Build

- Add gdal to devbox environment
- Update Devbox packages and environment variables

## [0.0.2] - 2025-12-31

### <!-- 0 -->⛰️  Features

- Introduce comprehensive build system and dependency management
- Convert geotiv to a header-only library

### <!-- 2 -->🚜 Refactor

- Migrate to new concord ecosystem (datapod, optinum, concord)

### <!-- 7 -->⚙️ Miscellaneous Tasks

- Clean up devbox and doctest configurations

### Build

- Add gdal to devbox environment
- Update Devbox packages and environment variables

## [3.1.0] - 2025-12-15

### Build

- Refactor build system to support both `xmake` and `cmake`

## [3.0.0] - 2025-11-01

### <!-- 0 -->⛰️  Features

- Refactor I/O into a library and introduce raster management
- Introduce 3D layer I/O and GeoTIFF parser utilities

### <!-- 7 -->⚙️ Miscellaneous Tasks

- Streamline build process and dependencies

## [3.0.0] - 2025-11-01

### <!-- 0 -->⛰️  Features

- Introduce 3D layer I/O and GeoTIFF parser utilities

### <!-- 7 -->⚙️ Miscellaneous Tasks

- Streamline build process and dependencies

## [2.1.0] - 2025-08-19

### <!-- 0 -->⛰️  Features

- Enhance GeoTIFF layer handling and examples

### Build

- Update build configuration for `concord` and CMake

## [2.0.2] - 2025-08-08

### <!-- 0 -->⛰️  Features

- Refactor: Improve coordinate system conversions

## [2.0.1] - 2025-08-08

### <!-- 1 -->🐛 Bug Fixes

- Adjust pixel data handling for accurate representation
- Adjust GeoTIFF image and coordinate system metadata

### <!-- 4 -->⚡ Performance

- Improve GeoTIFF pixel scale conversions with `concord`

## [2.0.0] - 2025-07-09

### <!-- 0 -->⛰️  Features

- Improve GeoTIFF CRS handling and WGS84 compatibility

## [1.2.0] - 2025-07-09

### <!-- 7 -->⚙️ Miscellaneous Tasks

- Configure CI for Docker image builds

## [1.1.4] - 2025-07-03

### <!-- 2 -->🚜 Refactor

- Refactor: Update dependencies and simplify data access

## [1.1.3] - 2025-06-29

### <!-- 0 -->⛰️  Features

- Pass custom tags to raster layers

## [1.1.2] - 2025-06-29

### <!-- 0 -->⛰️  Features

- Integrate global properties with TIFF ASCII tags

## [1.1.1] - 2025-06-28

### <!-- 7 -->⚙️ Miscellaneous Tasks

- Update concord to version 2.0.4 and remove datum parameter

## [1.1.0] - 2025-06-28

### <!-- 0 -->⛰️  Features

- Refactor Raster class for improved grid name retrieval
- Add Raster class for grid layer management

## [1.0.0] - 2025-06-15

### <!-- 0 -->⛰️  Features

- Refactor codebase to use internal CRS type

### <!-- 3 -->📚 Documentation

- Clean up blank lines in documentation

## [0.3.0] - 2025-06-11

### <!-- 0 -->⛰️  Features

- Add examples for tiff image generation
- Add support for per-IFD metadata and custom tags

## [0.2.1] - 2025-06-11

### <!-- 1 -->🐛 Bug Fixes

- Rename geotiff namespace to geotiv

## [0.2.0] - 2025-06-11

### <!-- 0 -->⛰️  Features

- Improve raster writing performance and examples
- Add example for writing a checkerboard to GeoTIFF
- Add support for floating-point TIFF resolution values
- Add GeoTIFF creation and writing functionality
- Add initial GeoTIFF parser implementation

### <!-- 1 -->🐛 Bug Fixes

- Update TIFF file header tests and seeking
- Correct typo in geotiff namespace

### <!-- 2 -->🚜 Refactor

- Enhance GeoTIFF parsing and validation for reading
- Add support for multi-layer GeoTIFFs
- Refactor TIFF reading and writing with Layer struct
- Refactor geotiff parsing and representation

### <!-- 3 -->📚 Documentation

- Update project and documentation
- Update logo files
- Improve documentation

### <!-- 4 -->⚡ Performance

- Improve GeoTIFF writing for raster collections

### <!-- 5 -->🎨 Styling

- Update logo color throughout the UI

### <!-- 7 -->⚙️ Miscellaneous Tasks

- Install changelog generator cliff

### Build

- Add new function run and update build aliases
- Reorganize geotiff parser and enable examples
- Initialize project with CMake build system
- Prepare development environment for direnv

## [1.3.7] - 2025-06-11

### <!-- 1 -->🐛 Bug Fixes

- Update TIFF file header tests and seeking

## [1.3.6] - 2025-06-11

### <!-- 7 -->⚙️ Miscellaneous Tasks

- Install changelog generator cliff

All notable changes to this project will be documented in this file.

