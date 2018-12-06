module MyModule

using LinearAlgebra

#Creates Meshgrid
include("meshgrid.jl")
include("normr.jl")

#Analytical functions (not sped these up)
include("Barber1992_GlideDislocation.jl")
include("Okada1985RectangularDislocation.jl")

#Continuum mechanics funcs
include("CalculateDSandSSDirs.jl")
include("RotateCosine3D.jl")
include("CheckDirectionCosinePerp.jl")
include("CreateDirectionCosinesFromStrikeAndDip.jl")

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

#Drawing 
include("contourfill.jl")

end