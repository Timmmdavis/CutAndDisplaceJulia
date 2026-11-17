#
# Correctness Tests
#

using Test 
using CutAndDisplaceJulia

fatalerrors = length(ARGS) > 0 && ARGS[1] == "-f"
quiet = length(ARGS) > 0 && ARGS[1] == "-q"
anyerrors = false


my_tests = ["test_TDvsBurgmannLinearFrictionCrack.jl",
            "test_TDvsOkada.jl",
            "test_LDvsBarberGlideDisc.jl",
			"test_ElasticConstantsCheck.jl",
			"test_LDvsOkada.jl",
			"test_TDvsEshelbyPennyCrack.jl",
            "test_TDvsMogi.jl",
            "test_TDvsSavageGravityValley.jl",
            "test_TDvsStressIntensity.jl",
            "test_TractionFreeSurface.jl",
            "test_TDvsStressIntensityElliptical.jl",
            "test_TDvsStressIntensityEllipticalGradient.jl",
            "test_TDvsStressIntensityGradientPenny.jl",
            # Unit tests of separate parts of the code. Figures, and the CGAL and plotting
            # tests, only when run through UnitTests/RunUnitTests.jl
            "UnitTests/ut_01_MeshReadWrite.jl",
            "UnitTests/ut_02_MeshGeometryAndTipElements.jl",
            "UnitTests/ut_03_EdgeTriCleanup.jl",
            "UnitTests/ut_04_RemeshingCGAL.jl",
            "UnitTests/ut_05_StressAndTractions.jl",
            "UnitTests/ut_06_FluidVolumeSolve.jl",
            "UnitTests/ut_07_StressIntensityAndPropagation.jl",
            "UnitTests/ut_08_PlottingCrackViews.jl",
            "UnitTests/ut_09_StallCheckAndExport.jl",
            "UnitTests/ut_10_DipAzimuthAndNeighbours.jl"]
			
println("Running tests:")

for my_test in my_tests
    try
        println("\t\033[1m\033[32mStarting\033[0m: $(my_test)")
        include(my_test)
        println("\t\033[1m\033[32mPASSED\033[0m: $(my_test)")
    catch e
        global anyerrors = true
        println("\t\033[1m\033[31mFAILED\033[0m: $(my_test)")
        if fatalerrors
            rethrow(e)
        elseif !quiet
            showerror(stdout, e, backtrace())
            println()
        end
    end
end

if anyerrors
    throw("Tests failed")
end