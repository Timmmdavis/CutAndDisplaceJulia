module CutAndDisplaceJulia
using LinearAlgebra
using Profile
using DelimitedFiles
using Statistics

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

#Types to contain boundary conditions
mutable struct Stresses;	σxx;σyy;σzz;σxy;σxz;σyz; end
mutable struct Strains;		εxx;εyy;εzz;εxy;εxz;εyz; end 
mutable struct Tractions;	Tn;Tss;Tds; 			 end 	
mutable struct MixedBoundaryConditions;Stresses;Tractions;  end
export
    Stresses,
	Strains,
	Tractions,
	MixedBoundaryConditions

#Init some types for influence matricies
mutable struct TractionInf; 	
	DnTn;	DnTss;	DnTds	
	DssTn;	DssTss;	DssTds
	DdsTn;	DdsTss;	DdsTds
	
end 	
mutable struct DispInf; 	
	DnUx;	DnUy;	DnUz	
	DssUx;	DssUy;	DssUz
	DdsUx;	DdsUy; 	DdsUz
end 

export
    StressInf,
	DispInf


#New Guys 07/02/2019
include("CalculateNormalTraction3D.jl")
include("TractionVectorCartesianComponents3D.jl")
include("CalculateTractionInChosenDirection3D.jl")
include("CreateTriangleNormal.jl")
include("HookesLaw3DStrain2Stress.jl")
include("CreateMidPoint.jl")
include("ElasticConstantsCheck.jl")
include("CreateFaceNormalAndMidPoint.jl")
include("CreateP1P2P3.jl")
include("CalculateNormalAndShearTractions3D.jl")
include("RemoveFixedElements3D.jl")


include("SlipCalculator3D.jl")
include("RepeatStruct.jl")
include("SetupTractionVector.jl")
include("LoadData.jl")
include("DataAppender3D.jl")
include("CalculateInfluenceMatrices3D.jl")
include("KnockOutFixedRows.jl")

#Surface loading functions
include("GoCadAsciiReader.jl")
include("STLReader.jl")

#Creates Meshgrid
include("meshgrid.jl")
include("normr.jl")

#Analytical functions (not sped these up)
include("Barber1992_GlideDislocation.jl")
include("Okada1985RectangularDislocation.jl")
include("Eshelby1957_PennyCrackSlipProfile.jl")
include("Mogi1958_SphericalCavity.jl")
include("Savage1984_GravityValleyStress.jl")

#Continuum mechanics funcs
include("CalculateSSandDSDirs.jl")
include("RotateCosine3D.jl")
include("CheckDirectionCosinePerp.jl")
include("CreateDirectionCosinesFromStrikeAndDip.jl")
include("HookesLaw2DStrain2Stress.jl")
include("HookesLaw2DStress2Strain.jl")
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

#Simple handy functions
include("BAsPercentOfA.jl")

#LinearAlgebra funcs
include("cross!.jl")

#Mirrors of some MATLAB functions
include("cart2pol.jl")

end