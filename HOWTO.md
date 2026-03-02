# ElasticSim — Installation Guide

## Overview

ElasticSim is a discrete elastic shell simulator built on
[geometry-central](https://geometry-central.net/) and
[Polyscope](https://polyscope.run/). It uses CMake as its build system.

---

## Prerequisites

| Dependency | Role | How it is provided |
|---|---|---|
| **geometry-central** | Surface mesh, intrinsic geometry, half-edge data structures | Git submodule (`deps/geometry-central`, branch `DDG` from `nzfeng/geometry-central`) |
| **Polyscope** | 3-D mesh visualization and GUI | Git submodule (`deps/polyscope`, branch `DDG` from `nzfeng/polyscope`) |
| **Google Test** | Unit testing framework | Git submodule (`deps/googletest`) |
| **Eigen 3.3+** | Linear algebra (vectors, matrices) | Auto-resolved by geometry-central (system install, or auto-downloaded) |
| **GNU Scientific Library (GSL)** | `gsl_multimin` conjugate-gradient optimizer | Must be installed separately (see below) |
| **CMake >= 3.10** | Build system | Must be installed separately |
| **C++17 compiler** | GCC >= 9, Clang >= 10, or MSVC >= 19.20 | Must be installed separately |
| **OpenGL + GLFW** | Required by Polyscope for rendering | System packages (see below) |

### Note on the submodule forks

Both `geometry-central` and `polyscope` point to the `DDG` branch of
`nzfeng`'s forks, not to the upstream repositories. These forks include
additional surface-mesh features (e.g. `DependentQuantityD`) used by the
elastic geometry classes. Using a different branch or the upstream repos will
cause build failures.

---

## Step 1 — Clone the repository and initialize submodules

```bash
git clone <repository-url> ElasticSim
cd ElasticSim
git submodule update --init --recursive
```

The `--recursive` flag is important because Polyscope itself has nested
submodules (GLFW, GLM, Dear ImGui, glad, stb, happly).

---

## Step 2 — Install external dependencies

### Linux (Debian / Ubuntu)

```bash
sudo apt-get update
sudo apt-get install build-essential cmake libgsl-dev \
    libgl-dev libxrandr-dev libxinerama-dev libxcursor-dev libxi-dev
```

The `libxrandr/xinerama/xcursor/xi` packages are required by GLFW (built as
part of Polyscope).

### macOS (Homebrew)

```bash
brew install cmake gsl
```

OpenGL and the windowing backend are provided by the OS frameworks.

### Windows (vcpkg)

```powershell
vcpkg install gsl:x64-windows
```

When running CMake, pass the vcpkg toolchain file:

```powershell
cmake .. -DCMAKE_TOOLCHAIN_FILE=<vcpkg-root>/scripts/buildsystems/vcpkg.cmake
```

---

## Step 3 — Required modifications to CMakeLists.txt

The project CMakeLists at `projects/elastic_membrane/CMakeLists.txt` ships
with **hardcoded absolute paths** that are specific to the original developer's
Windows machine. You **must** edit these before building.

### 3a. Fix the `add_subdirectory` calls for geometry-central and polyscope

**Original (lines 86–87):**
```cmake
add_subdirectory(C:/Users/dgrossma/Documents/GitHub/geometry-central deps/geometry-central)
add_subdirectory(C:/Users/dgrossma/Documents/GitHub/polyscope deps/polyscope)
```

**Change to** (use the repo's own submodules):
```cmake
add_subdirectory(../../deps/geometry-central deps/geometry-central)
add_subdirectory(../../deps/polyscope deps/polyscope)
```

### 3b. Fix the GSL include/library paths

**Original (lines 99–101):**
```cmake
target_include_directories(main PUBLIC "C:/Users/dgrossma/Documents/DEV/vcpkg/packages/gsl_x64-windows/include")
target_include_directories(main PUBLIC "C:/Users/dgrossma/Documents/DEV/vcpkg/packages/gsl_x64-windows/bin")
include("C:/Users/dgrossma/Documents/DEV/vcpkg/scripts/buildsystems/vcpkg.cmake")
```

**Change to** (let CMake's `find_package` handle it):
```cmake
# Delete the three lines above entirely.
# The existing find_package(GSL REQUIRED) on line 102 is sufficient
# when GSL is installed system-wide (Linux/macOS) or when the vcpkg
# toolchain file is passed via the cmake command line (Windows).
```

On Linux/macOS with GSL installed via the package manager, no extra include
paths are needed — `find_package(GSL REQUIRED)` does the right thing.

On Windows with vcpkg, pass `-DCMAKE_TOOLCHAIN_FILE=...` on the command line
instead of hardcoding it in CMakeLists.txt (see Step 2).

### 3c. Fix the args.hxx include path

**Original (line 104):**
```cmake
target_include_directories(main PUBLIC "C:/Users/dgrossma/Documents/GitHub/polyscope/deps/args")
```

**Change to:**
```cmake
target_include_directories(main PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args")
```

### 3d. (Optional) Update the C++ standard

The CMakeLists sets `-std=c++11`, but the source code uses `<format>` (C++20)
and `<Windows.h>` (Windows-only). If building on Linux/macOS:

- Replace or remove the `#include <Windows.h>` in `src/main.cpp` (it is only
  used for `Sleep()` — substitute with `std::this_thread::sleep_for` from
  `<thread>`).
- Either replace `std::format` usage with `printf`/`sprintf`, or bump the
  standard to C++20:

```cmake
SET(BASE_CXX_FLAGS "-std=c++20 -Wall -Wextra")
```

---

## Step 4 — Build

```bash
cd projects/elastic_membrane
mkdir build && cd build
cmake ..
make -j$(nproc)
```

The executable will be placed in `build/bin/main`.

On Windows with Visual Studio, open the folder in VS (it reads
`CMakeSettings.json`) or use:

```powershell
cd projects\elastic_membrane
mkdir build && cd build
cmake .. -G "Visual Studio 17 2022" -A x64 -DCMAKE_TOOLCHAIN_FILE=<vcpkg-root>/scripts/buildsystems/vcpkg.cmake
cmake --build . --config Release
```

---

## Step 5 — Run

The executable expects an input mesh (OBJ format). Several sample meshes are
provided in the `input/` directory at the repository root:

```bash
./bin/main ../../input/sphere.obj
```

---

## Summary of CMakeLists.txt changes

| Line(s) | What to change | Why |
|---|---|---|
| 86–87 | Point `add_subdirectory` to `../../deps/geometry-central` and `../../deps/polyscope` | Original uses hardcoded Windows paths |
| 99–101 | Delete the three hardcoded vcpkg lines | `find_package(GSL)` handles this portably |
| 104 | Use `${CMAKE_CURRENT_SOURCE_DIR}/../../deps/polyscope/deps/args` | Original uses hardcoded Windows path |
| 24 | (Optional) Change `-std=c++11` to `-std=c++20` | Source uses C++20 `<format>` |

---

## Project structure

```
ElasticSim/
├── deps/
│   ├── geometry-central/   # git submodule (nzfeng/geometry-central, DDG branch)
│   ├── polyscope/          # git submodule (nzfeng/polyscope, DDG branch)
│   └── googletest/         # git submodule
├── input/                  # sample OBJ meshes (sphere, torus, bunny, etc.)
├── utils/
│   ├── include/            # colormap.h, distortion.h, setup.h, solvers.h
│   └── src/                # colormap.cpp, distortion.cpp, solvers.cpp
└── projects/
    └── elastic_membrane/
        ├── CMakeLists.txt  # ← the file you need to edit
        └── src/
            ├── main.cpp
            ├── ElasticGeometry.h / .cpp
            └── ElasticGeomterySphericalCoor.h / .cpp
```
