module MyModule

using LinearAlgebra

#Creates Meshgrid
include("meshgrid.jl")
include("normr.jl")


#Continuum mechanics funcs
include("CalculateDSandSSDirs.jl")


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
include("TDdispHS.jl")
include("TDSetupD.jl")
include("AngDisDisp.jl")



#halfspace
include("TDdisp_HarFunc.jl")
include("AngSetupFSC.jl")
include("AngDisDispFSC.jl")

#Strain
include("TDstrainFS.jl")
include("TDSetupS.jl")
include("AngDisStrain.jl")
include("TensorTransformation3D.jl")
include("TensorTransformation2D.jl")

#Rotation in Mehdis funcs
include("RotateObject3DNewCoords.jl") #Mehdis 'CoordTrans'
include("RotateObject2D.jl")
#Good for both
include("trimodefinder.jl")
include("CalculateLocalTriCoords.jl")
include("CalcTDVectsAndAngles.jl")
include("TransformToADCS.jl")

include("Okada1985RectangularDislocation.jl")





end