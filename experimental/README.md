# experimental/ (not loaded by the package, not run by the tests)

Code and tests from 17-19 Nov 2026 that were taken out of the build so that it matches the Tuesday 17 Nov
08:55 solver (commit `a120c90`). Nothing here is deleted: the full history is on GitHub (`master` at
`d0e338c`, tag `backup-main-2026-11-19`).

| File | What it was |
|---|---|
| `FluidPressureSolver.jl` | Fluid volume solve inside the friction LCP (called from the rewritten `SlipCalculator3D`, `d30a264`) |
| `GradedRemesh.jl` | Graded mesh: small equilateral tip triangles on the front growing inwards (`79b39e6`) |
| `test/ut_03_EdgeTriCleanup.jl` | Unit test of the rewritten `IsosceliseEdgeTris` (`d30a264`, `651272e`) |
| `test/ut_06_FluidVolumeSolve.jl`, `test/ut_12_FluidSolverReduced.jl` | Unit tests of the new fluid volume solver |
| `test/ut_11_GradedRemesh.jl` | Unit test of the graded remesh |
| `test/test_TDvsStressIntensityEllipticalFixedVolume.jl` | Fixed volume stress intensity check using the new solver and graded mesh |

To try them again, move the source files back into `src/`, add them to `src/CutAndDisplaceJulia.jl`, and
restore the rewritten `SlipCalculator3D.jl` and `IsosceliseEdgeTris.jl` from commit `d0e338c`.
