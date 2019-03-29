#
# Correctness Tests
#

using Test 
using CutAndDisplaceJulia

fatalerrors = length(ARGS) > 0 && ARGS[1] == "-f"
quiet = length(ARGS) > 0 && ARGS[1] == "-q"
anyerrors = false

my_tests = ["test_TDvsOkada.jl",
            "test_LDvsBarberGlideDisc.jl",
			"test_ElasticConstantsCheck.jl",
			"test_LDvsOkada.jl",
			"test_TDvsEshelbyPennyCrack.jl",
            "test_TDvsMogi.jl",
            "test_TDvsSavageGravityValley.jl"]
			
			
println("Running tests:")

for my_test in my_tests
    try
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