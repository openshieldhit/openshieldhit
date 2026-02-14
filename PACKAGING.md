# OpenShieldHIT Packaging Guide

This document describes how to build and distribute OpenShieldHIT packages.

## Overview

OpenShieldHIT uses CMake and CPack for building distributable packages. The version number is automatically extracted from git tags and commit hashes.

### Supported Package Formats

- **Linux**: DEB (Debian/Ubuntu) and TGZ (source/binary tarball)
- **macOS**: TGZ (tarball)
- **Windows**: ZIP archive

## Version Information

The version is automatically generated at build time from git information:

```bash
git describe --tags --always --dirty
```

This produces version strings like:
- `v1.2.3` - Release version
- `v1.2.3-5-gabcd123` - 5 commits after v1.2.3
- `v1.2.3-dirty` - Uncommitted changes
- `abcd123` - In a non-tagged repository

To create a release version, create an annotated git tag:

```bash
git tag -a v1.0.0 -m "Release version 1.0.0"
git push origin v1.0.0
```

## Building Packages Locally

### Prerequisites

On **Ubuntu/Debian**:
```bash
sudo apt-get install cmake build-essential libsdl2-dev
```

On **macOS**:
```bash
brew install cmake sdl2
```

On **Windows**:
- Install CMake (https://cmake.org/download/)
- Install Visual Studio Build Tools or Visual Studio Community
- Install SDL2 development libraries

### Building a Release Package

#### 1. Configure and Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DOSH_BUILD_EXAMPLES=ON
cmake --build build
```

#### 2. Create DEB Package (Linux only)

```bash
cmake --build build --target package
```

This generates `openshieldhit-<version>-Linux.deb` in the build directory.

#### 3. Create Tarball Package (All platforms)

```bash
cmake --build build --target package
```

This generates `openshieldhit-<version>-Linux.tar.gz` (or `.tar.Z` on some systems).

### Installing a Local Package

#### From DEB package

```bash
sudo dpkg -i openshieldhit-<version>-Linux.deb

# Or with apt if placed in a repository
sudo apt install ./openshieldhit-<version>-Linux.deb
```

#### From tarball

```bash
tar -xzf openshieldhit-<version>-Linux.tar.gz
cd openshieldhit-<version>-Linux
# Follow standard installation method
```

### Verifying Installed Package

```bash
# Test main executable
osh --version

# Test examples (if installed)
bnct_sdl /usr/local/share/openshieldhit/examples/02_bnct/geo_cell.dat

# Verify installed files
dpkg -L openshieldhit  # For DEB packages
```

## GitHub Releases

### Automatic Release Artifacts

When you push a git tag, the CI/CD pipeline automatically:

1. Builds the project on Linux, macOS, and Windows
2. Creates DEB and TGZ packages on Linux
3. Attaches DEB packages as release artifacts to the GitHub release

### Manual GitHub Release Creation

```bash
# Create and push a release tag
git tag -a v1.0.0 -m "Release v1.0.0"
git push origin v1.0.0

# GitHub Actions will:
# - Build packages automatically
# - Create a GitHub Release
# - Attach DEB/TGZ artifacts
```

Then:
1. Go to https://github.com/grzanka/openshieldhit/releases
2. The latest release will have DEB packages attached
3. Download and install:
   ```bash
   wget https://github.com/grzanka/openshieldhit/releases/download/v1.0.0/openshieldhit-v1.0.0-Linux.deb
   sudo apt install ./openshieldhit-v1.0.0-Linux.deb
   ```

### Release Workflow

The GitHub Actions workflow includes:

1. **Build & Test** (all platforms): Compiles and runs tests
2. **Package** (Linux only, on success): Creates DEB and tarball packages
3. **Test Package** (Linux only, on packaging success): Installs and tests the DEB
4. **Release** (on git tags): Attaches packages to GitHub Release

## Package Contents

### DEB Package includes:

- **Executable**: `/usr/local/bin/osh`
- **Examples** (if OSH_BUILD_EXAMPLES=ON):
  - `/usr/local/bin/gemca_sdl_viewer`
  - `/usr/local/bin/bnct_sdl`
  - Example data: `/usr/local/share/openshieldhit/examples/`
- **Documentation**: `/usr/local/share/openshieldhit/README.md`
- **License**: `/usr/local/share/openshieldhit/LICENSE`

### Intermediate Build Artifacts

The following are created during build but not included in packages:

- `.a` (static library files) - Only for internal use
- Version header `build/version.h` - Auto-generated
- Object files, CMake cache, etc. - Standard CMake artifacts

## Advanced Packaging Topics

### Custom Installation Prefix

```bash
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/opt/openshieldhit
cmake --build build --target install
```

### Build Without Examples

```bash
cmake -S . -B build -DOSH_BUILD_EXAMPLES=OFF
cmake --build build --target package
```

### Create Source Package

Current setup creates binary packages. To create source packages:

```bash
# Use git archive instead
git archive --prefix openshieldhit-v1.0.0/ -o openshieldhit-v1.0.0-src.tar.gz v1.0.0
```

## Troubleshooting

### DEB not found after build

Check that you're on Linux:
```bash
uname -s  # Should be "Linux"
```

### Version shows as "0.0.1-unknown"

Git is not available or not in a git repository:
```bash
git describe --tags --always --dirty  # Should return version tag
```

### SDL2 not found when building from package

Install SDL2 development files:
```bash
sudo apt-get install libsdl2-dev
```

Or disable examples when building:
```bash
cmake -S . -B build -DOSH_BUILD_EXAMPLES=OFF
```

## See Also

- CMake Documentation: https://cmake.org/documentation/
- CPack Documentation: https://cmake.org/cmake/help/latest/module/CPack.html
- Debian Packaging Guide: https://wiki.debian.org/HowToPackageForDebian
