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
include("HookesLaw2dStrain2Stress.jl")
include("HookesLaw2dStress2Strain.jl")
include("TensorTransformation3D.jl")
include("TensorTransformation2D.jl")

#Normal line dislocation funcs
include("LD.jl")
include("LD_sum.jl")

#Normal tri dislocation funcs
include("TD.jl")
include("TD_sum.jl")

#Displacement
#include("TDdispFS.jl")
#include("TDdispHS.jl")
#include("TDSetupD.jl")
#include("AngDisDisp.jl")

#Strain
#include("TDstrainFS.jl")
#include("TDstrain_HarFunc.jl")
#include("TDSetupS.jl")
#include("AngDisStrain.jl")
#include("AngDisStrainFSC.jl")
#include("AngSetupStrainFSC.jl")

#halfspace
#include("TDdisp_HarFunc.jl")
#include("AngSetupDispFSC.jl")
#include("AngDisDispFSC.jl")

#Rotation in Mehdis funcs
include("RotateObject3DNewCoords.jl") #Mehdis 'CoordTrans'
include("RotateObject3DNewCoords!.jl") #Mehdis 'CoordTrans'
include("RotateObject2D.jl")
include("RotateObject2D!.jl")



#Good for both
#include("trimodefinder.jl")
#include("CalculateLocalTriCoords.jl")
#include("CalcTDVectsAndAngles.jl")
#include("TransformToADCS.jl")
#include("CalcSideVec.jl")
#include("CalcSlipVectorDiscCoords.jl")


#Drawing 
include("contourfill.jl")

#DELETE
include("TestFooAllocs.jl")


end