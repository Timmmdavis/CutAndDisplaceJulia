#Dependencies
using Test 
using LinearAlgebra
#Add MyModule
using MyModule

println("Starting Test 1")
##Test 1, Check the TDE's match Okada dislocations (displacement HS). Need to add strain!!!
include("test_TDvsOkada.jl")
println("Test 1 passed")

println("Starting Test 2")
##Test 2, Check the LD's match glide dislocation formulas
include("test_LDvsBarberGlideDisc.jl")
println("Test 2 passed")

println("Starting Test 3")
##Test 2, Check the LD's match half space surface displacement and strain of Okada
include("test_LDvsOkada.jl")
println("Test 3 passed")