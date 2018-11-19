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
include("RotateObject3DNewCoords.jl") #Mehdis 'CoordTrans'
include("RotateObject2D.jl")
include("trimodefinder.jl")
include("TDSetupD.jl")
include("AngDisDisp.jl")
include("CalculateLocalTriCoords.jl")


using LinearAlgebra: cross,norm

##Testing vectorised vs loops in julia 1.02
#include("DotLoopTest.jl")
#include("DotLoopTestVect.jl")
#include("DotLoopTest2.jl")

end