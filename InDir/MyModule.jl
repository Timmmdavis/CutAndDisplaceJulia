module MyModule

#Creates Meshgrid
include("meshgrid.jl")

#Normal line dislocation funcs
include("LDstressHS.jl")
include("LDstressFS.jl")
include("LDdispFS.jl")
include("LDdispHS.jl")

#More involved ones for create inf matricies. 
include("LDFSInfMat.jl")
include("LDHSInfMat.jl")

#Normal tri dislocation funcs
include("TDdispFS.jl")
#include("CoordTrans.jl")
#include("trimodefinder.jl")
#include("TDSetupD.jl")
#include("AngDisDisp.jl")

using LinearAlgebra: cross,norm

#include("cross.jl")


end