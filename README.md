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

## Prerequisites

Before installing, ensure you have:
- Administrator access to your computer
- A path with no spaces (OneDrive and similar services may cause issues)

## 1. Install Julia

1. Download and install Julia 1.11.3 from the official website
2. **IMPORTANT**: Enable the "Add Julia to PATH" option during installation
3. Open Julia as administrator (right-click and select "Run as administrator")

## 2. Install Required Packages
In Julia, enter the package manager by typing `]` and run:
```julia
] add https://github.com/Timmmdavis/FischerNewton.git
] add https://github.com/Timmmdavis/CutAndDisplaceJulia.git
] add https://github.com/Timmmdavis/CutAndDisplaceJuliaPlots.git
] add https://github.com/Timmmdavis/BuildCGAL.git
```

Paste the following into your Julia shell; we do this as 7zip in Julia moved directory and CMake needs it in a certain place

<details>
<summary>Click to expand 7zip file copy script</summary>

```julia
#!/usr/bin/env julia
using Printf
function copy_julia_files()
    # Get Julia installation directory
    julia_install_dir = Sys.BINDIR |> dirname
    
    println("Found Julia installation directory: ", julia_install_dir)
    
    # Source directory
    source_dir = joinpath(julia_install_dir, "libexec", "julia")
    
    if !isdir(source_dir)
        error("Source directory not found: $source_dir")
    end
    
    # Target directories
    target_dirs = [
        joinpath(julia_install_dir, "bin"),
        joinpath(julia_install_dir, "libexec")
    ]
    
    # Check if running with admin privileges
    is_admin = false
    try
        # This is a simple test to check for admin privileges
        # Try to write to a file in a typically protected directory
        test_file = joinpath(julia_install_dir, ".admin_test")
        open(test_file, "w") do io
            write(io, "test")
        end
        rm(test_file)
        is_admin = true
    catch
        is_admin = false
    end
    
    if !is_admin
        error("This script requires administrator privileges to copy files to system directories.")
    end
    
    # Find all files in source directory
    source_files = readdir(source_dir, join=true)
    
    # Copy files to target directories
    println("Copying files from: $source_dir")
    
    for target_dir in target_dirs
        println("Copying to: $target_dir")
        
        # Create target directory if it doesn't exist
        if !isdir(target_dir)
            println("Creating directory: $target_dir")
            mkdir(target_dir)
        end
        
        # Copy each file
        files_copied = 0
        for source_file in source_files
            file_name = basename(source_file)
            target_file = joinpath(target_dir, file_name)
            
            try
                cp(source_file, target_file, force=true)
                files_copied += 1
                @printf("  Copied: %s\n", file_name)
            catch e
                @printf("  Failed to copy %s: %s\n", file_name, e)
            end
        end
        
        println("  Successfully copied $files_copied files to $target_dir")
    end
    
    println("\nOperation completed!")
end
# Run the function
try
    copy_julia_files()
catch e
    println("Error: ", e)
    println("\nIf you're not running as administrator, please restart Julia with admin privileges and try again.")
end
```
</details>

Then 
```julia
] add https://github.com/Timmmdavis/BuildCGAL.git
#If that fails try
using Pkg
Pkg.build("CMake")
#Before repeating add BuildCGAL
```

## 3. Verify Your Installation

Check that everything is working correctly:

```julia
using FischerNewton
using CutAndDisplaceJulia
using CutAndDisplaceJuliaPlots
using BuildCGAL
```

If no errors appear, your installation was successful.

## 4. Test the Packages

Run the test suites to ensure all functionality works as expected:

```julia
] test FischerNewton
] test CutAndDisplaceJulia
```

For CutAndDisplaceJuliaPlots, use the test script from the GitHub readme.

## 5. Set Up CGAL Support

CGAL is essential for mesh creation and manipulation, particularly for dyke propagation simulations. If you are on Linux, this is going to be easy: see https://github.com/Timmmdavis/BuildCGAL (CGAL-4.13.1)

### Install 7-Zip
Test is 7-Zip is installed 
```julia
println(try output = read(`7z`, String); "7zip is installed: " * strip(match(r"7-Zip\s+[\d\.]+(\s+\([^\)]+\))?", output).match) catch; try output = read(Sys.iswindows() ? `"C:\\Program Files\\7-Zip\\7z.exe"` : `7za`, String); "7zip is installed: " * strip(match(r"7-Zip\s+[\d\.]+(\s+\([^\)]+\))?", output).match) catch; "7zip is NOT installed" end end)
```
If not
1. Download 7-Zip (64-bit, version 19.00+) from 7-zip.org
2. Add to your system PATH:
   - Open "Environment Variables" in Windows
   - Edit PATH and add `C:\Program Files\7-Zip\`
   - Restart Julia

### Configure Julia Environment

- Locate your Julia installation directory
- Copy files from `JULIA_INSTALL_DIR\libexec\julia` to both:
  - `JULIA_INSTALL_DIR\bin`
  - `JULIA_INSTALL_DIR\libexec`

## 6. Configure BuildCGAL

1. Install the Revise package:
```julia
] add Revise
using Revise
using BuildCGAL
```

2. Install dependencies:
```julia
include("C:/Users/USERNAME/.julia/packages/BuildCGAL/src/DownloadAndCheckAllCGALCompilerDependencies.jl")
```
   - You may need to run this multiple times
   - If links are broken, see the Troubleshooting section

3. Add Boost to your environment variables (as instructed by the installer)

### Compile CGAL Functions

```julia
BuildCGAL.Compile_CGAL_Meshing_Functions()
```

Test the compilation:
```julia
cd("desired_output_path")
] test BuildCGAL
```

Verify you see reasonable mesh outputs in your specified directory.

## 7. Running a Dyke Propagation Simulation

Navigate to your desired output directory and run:

```julia
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

## Troubleshooting

If you encounter issues with the CGAL installation, try these alternative approaches:

### Alternative 1: Manual CMake Installation

```julia
] add CMake
] precompile
] build CMake

# Add to PATH
# C:\Users\USERNAME\.julia\packages\CMake\ULbyn\deps\usr\bin

# Force source build if needed
ENV["CMAKE_JL_BUILD_FROM_SOURCE"] = 1
] build CMake
```

### Alternative 2: Manual MSYS2/MinGW Installation

1. Download MSYS2 from msys2.org
2. Open MSYS2 MINGW64 terminal and run:
```
pacman -S mingw-w64-i686-toolchain
```
3. Verify with:
```
gcc --version
```

### CMake Generator Issues

If the MinGW Makefiles generator isn't found:
```
cmake --help
```
Try alternative generators like "MSYS Makefiles" or "NMake Makefiles"

### Package Version Conflicts

Create a local environment:
```julia
] activate .
] add ./FischerNewton
] add ./CutAndDisplaceJulia
] add ./CutAndDisplaceJuliaPlots
] add ./BuildCGAL
```

Or use development mode:
```julia
] dev ./FischerNewton
] dev ./CutAndDisplaceJulia
] dev ./CutAndDisplaceJuliaPlots
] dev ./BuildCGAL
```

You may need to adjust Project.toml files manually for compatibility.

### General Troubleshooting Tips

- Ensure paths have no spaces and are properly quoted
- Verify all environment variables are correctly set
- Check that dependencies are in your PATH
- For specific errors, search online using the exact error message
