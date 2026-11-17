# CutAndDisplaceJulia

[![Build Status](https://travis-ci.com/Timmmdavis/CutAndDisplaceJulia.svg?token=1HhESyMNyqzV8R22Pqq6&branch=master)](https://travis-ci.com/Timmmdavis/CutAndDisplaceJulia)
[![Coverage Status](https://codecov.io/gh/Timmmdavis/CutAndDisplaceJulia/branch/master/graph/badge.svg?token=IbbZ8n4385)](https://codecov.io/gh/Timmmdavis/CutAndDisplaceJulia)
[![DOI](https://zenodo.org/badge/155747027.svg)](https://zenodo.org/badge/latestdoi/155747027)

![FractureAscendingDueToStressGradients](https://github.com/Timmmdavis/HostingGifs/blob/master/MeshPropBothAnglesFast.gif)  

## Overview

CutAndDisplaceJulia is a Julia implementation of the [CutAndDisplace](https://github.com/Timmmdavis/CutAndDisplace) project. This package provides 3D linear elastic boundary element code for modeling geophysical problems. While the 2D functionality has not been implemented (except for analytical components), the 3D functionality is fully operational and highly efficient.

## Dependencies

The package relies on two key dependencies:
- **[CGAL](https://github.com/CGAL/cgal)**: For mesh generation and manipulation, connected through the [BuildCGAL](https://github.com/Timmmdavis/BuildCGAL) interface
- **[FischerNewton](https://github.com/Timmmdavis/FischerNewton)**: For handling non-linear solving of complementarity conditions

## Getting Started

### Running a Test Case

To get started with this package:

1. Open Julia with necessary packages (like Revise) imported
2. Run the following code (Windows example):
```julia
using CutAndDisplaceJulia
include(raw"C:\Users\UserName\...\CutAndDisplaceJulia\test\test_TDvsMogi.jl")
```

This approach is similar to the 'MainFrame.m' file in the MATLAB version.

### Finding Your Installation Path

To locate your installation path:
```julia
println(pathof(CutAndDisplaceJulia))
```

## Core Functionality

The main computational functions are located in the `src` directory:

- **LD.jl**: Computes stresses and displacements using analytical line dislocations
- **TD.jl**: Computes stresses and displacements using triangular dislocations

While TD.jl has been tested for type stability, LD.jl would benefit from further testing. Potential areas for improvement in TD.jl include memory management in the `CalculateLocalTriCoords` and `GlobalToTDECoords` functions, as well as the sections in 'FS' functions that are activated when the 'ImageFlag' is set to 1.

> **Note**: When debugging or modifying this code, we recommend removing `Threads.@threads` from the outer loop of these functions.

## Performance Benchmarks

### Halfspace Comparison (Triangular Dislocations)

Performance comparison between Julia and MATLAB implementations for commit 25dffe8 using functions from [TDstrainHS.m](https://github.com/Timmmdavis/CutAndDisplace/blob/master/3dCode/TDFunctions/TDstrainHS.m) (commit b07067a):

Test environment: Windows, Intel(R) Xeon(R) CPU X5472 @ 3.00GHz, 8 logical processors

| No of Tris | No of obs Points | MATLAB (seconds) | Julia (seconds) | Speedup |
|------------|:----------------:|:----------------:|:---------------:|:-------:|
| 2          | 50×50            | 0.37             | 0.06            | 5.9×    |
| 10         | 50×50            | 1.69             | 0.12            | 13.4×   |
| 100        | 50×50            | 17.12            | 0.87            | 19.6×   |
| 2500       | 50×50            | 425.94           | 21.38           | 19.9×   |
| 100        | 2×2              | 0.72             | 0.003           | 240×    |
| 100        | 10×10            | 4.09             | 0.03            | 120×    |
| 100        | 50×50            | 17.12            | 0.87            | 19.6×   |
| 100        | 100×100          | 59.84            | 3.51            | 17×     |

### Fullspace Comparison (Triangular Dislocations)

Performance comparison for fullspace calculations:

| No of Tris | No of obs Points | MATLAB (seconds) | Julia (seconds) | Speedup |
|------------|:----------------:|:----------------:|:---------------:|:-------:|
| 2          | 50×50            | 0.013            | 0.005           | 2.6×    |
| 10         | 50×50            | 0.072            | 0.013           | 5.5×    |
| 100        | 50×50            | 0.907            | 0.104           | 8.7×    |
| 2500       | 50×50            | 186.39           | 4.12            | 45×     |
| 100        | 2×2              | 0.057            | 0.005           | 11.4×   |
| 100        | 10×10            | 0.086            | 0.04            | 2.15×   |
| 100        | 50×50            | 0.907            | 0.104           | 8.7×    |
| 100        | 100×100          | 3.155            | 0.413           | 7.6×    |
| 4          | 2×2              | 0.003            | 0.0001307       | 23×     |
| 16         | 4×4              | 0.013            | 0.00049         | 27×     |
| 100        | 10×10            | 0.09             | 0.00722         | 12×     |
| 400        | 20×20            | 1.30             | 0.1031          | 13×     |
| 900        | 30×30            | 10.44            | 0.518           | 20×     |
| 2500       | 50×50            | 197.95           | 4.12            | 48×     |
| 3600       | 60×60            | 723.90           | 8.64            | 84×     |
| 6400       | 80×80            | 2906.40          | 27.58           | 105×    |

> **Note**: TD.jl produces all 3 influence matrices, effectively doing 3 times more work than reported in the benchmarks.

# Installation Guide

## 1. Install Julia

Install Julia 1.11 or newer, ideally with [juliaup](https://julialang.org/install/). You don't need admin rights, and your folders can be anywhere (OneDrive and paths with spaces are fine).

## 2. Install the packages

Start Julia, press `]` to enter the package manager, and run:

```julia
] add https://github.com/Timmmdavis/FischerNewton.git
] add https://github.com/Timmmdavis/CutAndDisplaceJulia.git
] add https://github.com/Timmmdavis/CutAndDisplaceJuliaPlots.git
] add https://github.com/Timmmdavis/BuildCGAL.git
```

Keep this order: FischerNewton is not in Julia's General registry, so it has to be added before CutAndDisplaceJulia.

This installs everything, including CGAL. BuildCGAL downloads prebuilt CGAL meshing tools for your machine (Windows x86_64, Linux x86_64 or Apple Silicon Mac) and shows a progress bar while it does. Nothing is compiled, and no system settings or environment variables are changed.

Prefer one line? This does the same thing:

```julia
import Pkg; Pkg.add([Pkg.PackageSpec(url = "https://github.com/Timmmdavis/$p.git") for p in ("FischerNewton", "CutAndDisplaceJulia", "CutAndDisplaceJuliaPlots", "BuildCGAL")])
```

## 3. Check it works

```julia
using FischerNewton, CutAndDisplaceJulia, CutAndDisplaceJuliaPlots, BuildCGAL
BuildCGAL.tools_dir()       # where the CGAL tools were installed
] test BuildCGAL            # runs every CGAL tool on sample meshes (a few seconds)
] test FischerNewton
] test CutAndDisplaceJulia
```

For CutAndDisplaceJuliaPlots, use the test script from its GitHub README.

## 4. Running a Dyke Propagation Simulation

Start a fresh Julia session, go to the folder where you want the results, and run:

```julia
# Results and figures are written here
cd("path/to/output/folder")   

#Then load packages we will use
using CutAndDisplaceJulia, CutAndDisplaceJuliaPlots, BuildCGAL

ModuleDir = dirname(dirname(pathof(CutAndDisplaceJulia)))
TestScript = joinpath(ModuleDir, "test", "test_PropagationTest.jl")
include(TestScript)

# Configuration
HSFlag = 0  # Constant
printstyled("Hs off \n", color=:cyan) 

# Elastic constants
G = ShearModulus(2.0e9)
ν = PoissonsRatio(0.25)
(K, E, λ, ν, G) = CutAndDisplaceJulia.ElasticConstantsCheck(G, ν)

# Material properties
ρrock = 2900
ρfluid = 2600
g = 9.81
Δρ = ((ρrock-ρfluid)*g)
Kc = 5e7  # 50 MPa √m

# Crack volume calculation
CrackVolume = real(((1-ν)/(16*G))*((9*π^4*Kc^8)/(complex(Δρ)^5))^(1/3))
NoTris = 300

# Run simulation
(PropFlag, maxX, minX, maxY, minY, maxZ, minZ) = testProp(HSFlag, ν, G, Δρ, Kc, CrackVolume, NoTris)
```

Results will be saved to your current directory. View the meshes using meshlab.exe or the included plotting tools.


## Docker

To get a ready-to-run Linux environment instead, use the image in [`docker/`](docker):

```bash
docker build -t cutanddisplace-julia docker
docker run -it --rm -v /path/to/local/folder:/data cutanddisplace-julia
```

## Troubleshooting

- **`BuildCGAL has no prebuilt CGAL tools for this platform`**: you're on a platform without prebuilt binaries (for example an Intel Mac). Build the tools yourself from [BuildCGAL's `cgal_tools/` folder](https://github.com/Timmmdavis/BuildCGAL/tree/master/cgal_tools), then set `ENV["BUILDCGAL_TOOLS_DIR"]` to the folder containing them.
- **`The prebuilt CGAL tools have not been downloaded yet`**: the download was interrupted. Run `] instantiate` or `] add https://github.com/Timmmdavis/BuildCGAL.git` again.
- **Upgrading from an old install**: earlier versions of BuildCGAL added Make, MinGW, Boost and CGAL folders to your PATH and needed a 7-Zip copy step. None of that is needed now; see the [BuildCGAL README](https://github.com/Timmmdavis/BuildCGAL#upgrading-from-the-old-installer).
- **Package version conflicts**: install into a separate project environment so other packages can't clash:
  ```julia
  ] activate CutAndDisplace   # creates ./CutAndDisplace/Project.toml
  ] add https://github.com/Timmmdavis/FischerNewton.git
  # ... then the other three as above
  ```
- **Working on the code itself**: use `] dev` with the local clones instead of `] add`.
