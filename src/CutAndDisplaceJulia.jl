module CutAndDisplaceJulia
using LinearAlgebra
using Profile

#New Guys 07/02/2019
include("CalculateNormalTraction3d.jl")
include("TractionVectorCartesianComponents3d.jl")
include("CalculateTractionInChosenDirection3d.jl")
include("HookesLaw3dStrain2Stress.jl")
include("CreateMidPoint.jl")


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

include("TensorTransformation3D!.jl")

#Normal line dislocation funcs
include("LD.jl")
include("LD_sum.jl")

#Normal tri dislocation funcs
include("TD.jl")
include("TD_sum.jl")

#Rotation in Mehdis funcs
include("RotateObject3DNewCoords.jl") #Mehdis 'CoordTrans'
include("RotateObject3DNewCoords!.jl") #Mehdis 'CoordTrans'
include("RotateObject2D.jl")
include("RotateObject2D!.jl")

#Drawing 
include("contourfill.jl")

#LinearAlgebra funcs
include("cross!.jl")

include("tests.jl")
end