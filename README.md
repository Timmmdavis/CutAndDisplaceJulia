# CutAndDisplaceJulia

[![Coverage Status](https://codecov.io/github/Timmmdavis/CutAndDisplaceJulia/graph/badge.svg?token=IbbZ8n4385)](https://codecov.io/github/Timmmdavis/CutAndDisplaceJulia)
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

`] test CutAndDisplaceJulia` skips the tests that need CGAL. To run all the unit tests, including the CGAL
meshing ones, and save their figures (in `test/UnitTests/figures`):

```julia
include(joinpath(pkgdir(CutAndDisplaceJulia), "test", "UnitTests", "RunUnitTests.jl"))
```

## 4. Running a Dyke Propagation Simulation in a Valley

A buoyant, fluid-filled crack is started below an incised valley and propagated until it reaches the valley
surface or is trapped. The rock is loaded with the plane-strain Savage et al. (1984) gravity solution for a
valley (no tectonic stress), and the fluid as tractions on the crack. Axes: x across the valley, y along it,
z up (z = 0 is the original flat surface, so z < 0 below it).

Start a fresh Julia session, go to the folder where you want the results, and run:

```julia
# Results and figures are written here
cd("path/to/output/folder")

# Load the packages and the test
using CutAndDisplaceJulia, CutAndDisplaceJuliaPlots, BuildCGAL
ModuleDir = dirname(dirname(pathof(CutAndDisplaceJulia)))
include(joinpath(ModuleDir, "test", "test_PropagationValleyStress.jl"))

# ---------------------------------------------------------------------------------------------
# Solver settings
# ---------------------------------------------------------------------------------------------
HSFlag = 0                 # 0 = full space, 1 = half space (free-surface correction in the BEM)
NoTris = 350               # approximate number of triangles in the crack mesh
Meshing = :graded          # :uniform = one edge length everywhere; :graded = small triangles along the crack front,
                           #   growing inwards, for about the same number of triangles

# ---------------------------------------------------------------------------------------------
# Physical parameters of the rock and the fluid
# ---------------------------------------------------------------------------------------------
G = ShearModulus(2.0e9)    # shear modulus of the rock [Pa]
ν = PoissonsRatio(0.25)    # Poisson's ratio of the rock, used by the BEM (the valley loading has its own, see νLoad)
(K, E, λ, ν, G) = CutAndDisplaceJulia.ElasticConstantsCheck(G, ν)
Kc = 5e7                   # fracture toughness of the rock [Pa√m] (50 MPa√m); the crack grows where KI reaches this

# ---------------------------------------------------------------------------------------------
# Body forces (gravity and the densities that set the buoyancy)
# ---------------------------------------------------------------------------------------------
ρrock = 2900               # rock density [kg/m³]; gives the lithostatic stress ρrock*g*|z|
ρfluid = 2600              # fluid (magma) density [kg/m³]; its weight acts on the crack walls
g = 9.81                   # gravitational acceleration [m/s²]
Δρ = (ρrock - ρfluid)*g    # buoyancy gradient [Pa/m]; sets the critical radius and critical volume

# ---------------------------------------------------------------------------------------------
# The valley (Savage et al. 1984): surface z = a²b/(u²+a²) in the conformal coordinate u
# ---------------------------------------------------------------------------------------------
ValleyWidth = 2000.0       # a [m]: width of the valley (larger = wider, gentler valley)
ValleyDepth = -1000.0      # b [m]: depth of the valley floor below the original surface, NEGATIVE. Needs b > -a
νLoad = 0.5                # Poisson's ratio of the stress solution only (0.5 = isotropic far field, σxx = σzz).
                           #   The out-of-plane stress is σyy = νLoad*(σxx + σzz) (plane strain)

# ---------------------------------------------------------------------------------------------
# Where the crack starts (all distances in units of the valley width a)
# ---------------------------------------------------------------------------------------------
XFactor = 1.0              # crack centre at x = XFactor*a: horizontal distance from the valley axis (x = 0)
DepthFactor = 0.5          # crack centre at z = b - DepthFactor*a: distance below the valley floor level
# The crack starts as a penny whose plane contains the most compressive principal stress at its centre
# (its normal is the least compressive direction), so it opens against the least compressive stress.

# ---------------------------------------------------------------------------------------------
# How big the crack is and how much fluid it holds (relative to the flat, pre-valley critical values:
# a vertical crack under the lithostatic gradient, ρrock*g against ρfluid*g)
# ---------------------------------------------------------------------------------------------
RadiusFactor = 0.5         # start radius = RadiusFactor * critical radius (below 1: smaller than critical)
VolumeFactor = 1.2         # fluid volume = VolumeFactor * critical volume (below 1 the crack is likely trapped)

# ---------------------------------------------------------------------------------------------
# Run: stops when the crack reaches the valley surface (PropFlag = 1) or is trapped (PropFlag = 0)
# ---------------------------------------------------------------------------------------------
(PropFlag, maxX, minX, maxY, minY, maxZ, minZ) = testPropagationValleyStress(HSFlag, ν, G, Δρ, Kc, NoTris;
    ρrock = ρrock, ValleyWidth = ValleyWidth, ValleyDepth = ValleyDepth, νLoad = νLoad,
    XFactor = XFactor, DepthFactor = DepthFactor,
    RadiusFactor = RadiusFactor, VolumeFactor = VolumeFactor, Meshing = Meshing)
```

Before the first step a figure of the setup, `ValleyStress-Setup.png`, is written to the current folder. It
is in units of `a` and shows the valley surface, the most compressive stress (colour and bars) and the crack,
with a zoom on the crack so it is visible whatever its size.

`PropFlag` is 1 if the crack reached the valley surface, 0 if it was trapped (stopped growing), 2 if no part of
the front could propagate and 4 if the loop limit was reached. The crack is remeshed every step (see `Meshing`).
Results will be saved to your current directory. View the meshes using meshlab.exe or the included
plotting tools.

For a crack in a flat half space (a chosen dip and horizontal/vertical stress ratio instead of a valley), see
`testProp` in `test/test_PropagationTest.jl`.


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
