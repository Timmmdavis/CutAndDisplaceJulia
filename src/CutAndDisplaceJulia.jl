module CutAndDisplaceJulia
using LinearAlgebra
using Profile
using DelimitedFiles

#Init some types for elastics (ElasticConstantsCheck)
struct PoissonsRatio; 	ν::Float64;		end 
struct ShearModulus; 	G::Float64;		end
struct BulkModulus; 	K::Float64;		end
struct YoungsModulus; 	E::Float64;		end
struct LamesConstant; 	λ::Float64;		end
export
    ShearModulus,
	PoissonsRatio,
	BulkModulus,
	YoungsModulus,
	LamesConstant
	
#Init some types for influence matricies
struct StressInf; 	
	DnTn;	DnTss;	DnTds	
	DssTn;	DssTss;	DssTds
	DdsTn;	DdsTss;	DdsTds
	
end 	
struct DispInf; 	
	DnUx;	DnUy;	DnUz	
	DssUx;	DssUy;	DssUz
	DdsUx;	DdsUy; 	DdsUz
end 
export
    StressInf,
	DispInf


#New Guys 07/02/2019
include("CalculateNormalTraction3d.jl")
include("TractionVectorCartesianComponents3d.jl")
include("CalculateTractionInChosenDirection3d.jl")
include("CreateTriangleNormal.jl")
include("HookesLaw3dStrain2Stress.jl")
include("CreateMidPoint.jl")
include("ElasticConstantsCheck.jl")
include("CreateFaceNormalAndMidPoint.jl")
include("CreateP1P2P3.jl")

#Surface loading functions
include("GoCadAsciiReader.jl")
include("STLReader.jl")

#Creates Meshgrid
include("meshgrid.jl")
include("normr.jl")

#Analytical functions (not sped these up)
include("Barber1992_GlideDislocation.jl")
include("Okada1985RectangularDislocation.jl")

#Continuum mechanics funcs
include("CalculateSSandDSDirs.jl")
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