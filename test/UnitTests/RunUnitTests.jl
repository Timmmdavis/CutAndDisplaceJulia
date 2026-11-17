# Runs the unit tests in this folder AND draws their figures (test/UnitTests/figures/*.png).
#
# Needs CutAndDisplaceJulia, CutAndDisplaceJuliaPlots (for CairoMakie) and, for the remeshing
# tests, BuildCGAL in the active environment. Usage, from any folder:
#	julia RunUnitTests.jl            (all tests)
#	julia RunUnitTests.jl 07         (only tests whose file name contains "07")
#
# Inside Pkg.test (runtests.jl) the same files run without figures and without the BuildCGAL tests.

using CutAndDisplaceJulia
import CutAndDisplaceJuliaPlots       #(also the plotting under test, unit test 8)
const CM = CutAndDisplaceJuliaPlots.CairoMakie       #figures on (see UnitTestHelpers.jl)
try
	@eval import BuildCGAL                            #remeshing tests only run if this loads
catch
	println("BuildCGAL not available: the remeshing tests will be skipped")
end

using Test
Filter=length(ARGS)>0 ? ARGS[1] : ""
Files=sort(filter(f->startswith(f,"ut_") && endswith(f,".jl") && occursin(Filter,f),readdir(@__DIR__)))

@testset "CutAndDisplaceJulia unit tests" begin
	for f in Files
		println("\n== $f")
		t=@elapsed include(joinpath(@__DIR__,f))
		println("   ($(round(t,digits=1)) s)")
	end
end
