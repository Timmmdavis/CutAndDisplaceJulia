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
#Displacement
include("TDdispFS.jl")
include("TDSetupD.jl")
include("AngDisDisp.jl")
#Strain
include("TDstrainFS.jl")
include("TDSetupS.jl")
include("AngDisStrain.jl")
include("TensorTransformation3D.jl")
#Rotation in Mehdis funcs
include("RotateObject3DNewCoords.jl") #Mehdis 'CoordTrans'
include("RotateObject2D.jl")
#Good for both
include("trimodefinder.jl")
include("CalculateLocalTriCoords.jl")
include("CalcTDVectsAndAngles.jl")
include("TransformToADCS.jl")


using LinearAlgebra: cross,norm


end