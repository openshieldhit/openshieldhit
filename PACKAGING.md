# OpenShieldHIT Packaging Guide

## Version Management

Version is extracted from git at **CMake configure time** using `git describe --tags --always --dirty` and written into a generated header:

```
<build_dir>/generated/common/osh_version.h
```

This file is produced by `configure_file` from `src/common/osh_version.h.in` and provides the macros `OSH_VERSION`, `OSH_VERSION_MAJOR`, `OSH_VERSION_MINOR`, and `OSH_VERSION_PATCH` to the implementation. The public API exposes version information through functions (`openshieldhit_version_string()` etc.) rather than macros.

To update version in an existing build, reconfigure:
```bash
cmake -B build
```

## Creating Releases

1. **Create release manually on GitHub web**:
   - Go to https://github.com/openshieldhit/openshieldhit/releases
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

Install dependencies:
```bash
sudo apt-get install cmake build-essential
```

Build and package:
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
```

```bash
cmake --build build --target package
```

Install DEB package:
```bash
sudo apt-get install ./build/openshieldhit-*.deb
```

### macOS

Install dependencies:
```bash
brew install cmake
```

Build and package:
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
```

```bash
cmake --build build --target package
```

### Windows

Configure:
```bash
cmake -S . -B build
```

Build:
```bash
cmake --build build --config Release --target package
```

## Package Contents

**DEB package** (`/usr/local/`):
- `bin/openshieldhit` - Main executable
- `share/openshieldhit/README.md`, `LICENSE`

**Formats**: DEB (Linux), TGZ (all), ZIP (Windows)

## Verify Installation

Print version:
```bash
openshieldhit --version
```

Show help:
```bash
openshieldhit --help
```

List installed files (DEB packages):
```bash
dpkg -L openshieldhit
```

