# OpenShieldHIT Packaging Guide

## Version Management

Version is extracted from git at **CMake configure time** using `git describe --tags --always --dirty` and passed to the compiler as a preprocessor definition.

To update version in an existing build, reconfigure:
```bash
cmake -B build  # Re-runs configure, updates version
```

## Creating Releases

1. **Create release manually on GitHub web**:
   - Go to https://github.com/grzanka/openshieldhit/releases
   - Click "Create a new release"
   - Tag version: `v0.1.0`
   - Release title: `Release v0.1.0`
   - Click "Publish release"

2. **Workflow automatically starts**:
   - The **Package and Release** workflow detects the release publication
   - Builds DEB and TGZ packages on Linux
   - Tests the DEB package installation
   - Uploads both artifacts to your release

**Version format**: Use semantic versioning `vMAJOR.MINOR.PATCH` (e.g., `v0.1.0`, `v1.2.3`)

## GitHub Actions

### Automatic Triggers
- **Push to main branch**: Builds packages, uploads as artifacts (retention: 3 days)
- **Release published on GitHub web**: Builds packages and automatically uploads artifacts to the release

### Manual Trigger
1. Go to: https://github.com/openshieldhit/openshieldhit/actions
2. Select "Package and Release" workflow
3. Click "Run workflow" → Select branch → "Run workflow"

## Local Packaging

### Ubuntu/Debian
```bash
# Install dependencies
sudo apt-get install cmake build-essential libsdl2-dev

# Build and package
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DOSH_BUILD_EXAMPLES=ON
cmake --build build --target package

# Install
sudo dpkg -i build/openshieldhit-*.deb
```

### macOS
```bash
brew install cmake sdl2
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DOSH_BUILD_EXAMPLES=ON
cmake --build build --target package
```

### Windows
```bash
cmake -S . -B build -DOSH_BUILD_EXAMPLES=OFF
cmake --build build --config Release --target package
```

## Package Contents

**DEB package** (`/usr/local/`):
- `bin/osh` - Main executable
- `bin/gemca_sdl_viewer`, `bin/bnct_sdl` - Examples (if enabled)
- `share/openshieldhit/examples/` - Example data
- `share/openshieldhit/README.md`, `LICENSE`

**Formats**: DEB (Linux), TGZ (all), ZIP (Windows)

## Verify Installation

```bash
osh --version
osh --help
dpkg -L openshieldhit  # List package files
```

