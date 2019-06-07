module CutAndDisplaceJulia
using LinearAlgebra
using Profile
using DelimitedFiles
using Statistics
using FischerNewton
using Interpolations
using Random
using Statistics
using UnicodePlots
using DataFrames
using Printf

######################################################
println("Remove me!")
using Debugger
#=
#Assuming the func has a @bp statement
using Debugger
include(raw("functionpath"))
@enter function
C
c 'enter' #continue to bp
` #go into repl mode
ctrl+c #to exit back to debug
? #for help
=#
######################################################

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
struct Stresses;	σxx;σyy;σzz;σxy;σxz;σyz; end
struct Strains;		εxx;εyy;εzz;εxy;εxz;εyz; end 
struct Disps;		Ux;	Uy;	Uz; end 
struct Tractions;	Tn;Tss;Tds; 			 end 	
struct FrictionParameters; 	   µ;Sf;  		     	end
struct MixedBoundaryConditions;Stresses;Tractions;  end
struct MixedBoundaryConditionsFriction;MixedBoundaryConditions;FrictionParameters;	end
struct MixedBoundaryConditionsFluidVolume;MixedBoundaryConditions;Volumes;	end


println("I would prefer not be mutable")
mutable struct TriangleEdges;FeLe;FeMd;FeEv;FeM2Ev;FreeFlg;FeM2ELe;IntAng;K1;K2;K3;StrainEnergy; 	end

export
    Stresses,
	Strains,
	Disps,
	Tractions,
	FrictionParameters,
	MixedBoundaryConditions,
	MixedBoundaryConditionsFriction,
	MixedBoundaryConditionsFluidVolume,
	TriangleEdges

println("I would also prefer not be mutable")
#Init some types for influence matricies
mutable struct TractionInf; 	
	DnTn;	DnTss;	DnTds	
	DssTn;	DssTss;	DssTds
	DdsTn;	DdsTss;	DdsTds
end 	

struct StrainInf; 	
	εxxDn;	εyyDn; 	εzzDn;	εxyDn; 	εxzDn;	εyzDn 
	εxxDss;	εyyDss; εzzDss;	εxyDss; εxzDss;	εyzDss 
	εxxDds;	εyyDds; εzzDds;	εxyDds; εxzDds;	εyzDds 
end 	

struct DispInf; 	
	DnUx;	DnUy;	DnUz	
	DssUx;	DssUy;	DssUz
	DdsUx;	DdsUy; 	DdsUz
end 

mutable struct InfMat; A; end
mutable struct BoundaryConditionsVec; b; end

export
    TractionInf,
	DispInf,
	StrainInf,
	InfMat,
	BoundaryConditionsVec


include("OFFReader.jl")
include("xyzExport.jl")
include("OFFExport.jl")
include("ConnectedConstraints.jl")
include("qGetRotQuaternion.jl")
include("qRotatePoint.jl")
include("FindConnectedTriangles.jl")
include("ismember.jl")
include("MatchingRow.jl")
include("CreateSortedEdgePoints.jl")
include("IsosceliseEdgeTris.jl")
include("CleanEdgeTris.jl")
include("CreateTrianglesPointsFromP1P2P3.jl")
include("CalculateInternalTriAngles.jl")
include("RemoveDodgyNewEdges.jl")
include("ConnectedComponentsReader.jl")
include("GetSortedEdgesOfMeshList.jl")  
include("CollapseEdgeTris.jl") 
include("AngleBetweenVectors!.jl")  
include("FlipValue.jl")
include("FindIntersectionOf3DVectors.jl")
include("IsosceliseEdgeTrisNew.jl")
include("ExportCrackMesh.jl")
include("STLExport.jl")


#New Guys 07/02/2019
include("CalculateNormalTraction3D.jl")
include("CalculateNormalTraction3D!.jl")
include("TractionVectorCartesianComponents3D.jl")
include("TractionVectorCartesianComponents3D!.jl")
include("CalculateTractionInChosenDirection3D.jl")
include("CalculateTractionInChosenDirection3D!.jl")
include("CreateTriangleNormal.jl")
include("HookesLaw3DStrain2Stress.jl")
include("HookesLaw3DStrain2Stress!.jl")
include("CreateMidPoint.jl")
include("ElasticConstantsCheck.jl")
include("CreateFaceNormalAndMidPoint.jl")
include("CreateP1P2P3.jl")
include("CalculateNormalAndShearTractions3D.jl")
include("RemoveFixedElements3D.jl")
include("AreaOfTriangle3D.jl")
include("GetCrackTipElements3D.jl")
include("EdgeConstraints.jl")
include("RotateObject3DAllignVectors.jl")
include("StressIntensity3D.jl")
include("CalculateNormalShearTraction2D.jl")
include("Tada_StrIntInclinedPennyTension.jl")


#Fluid filled crack bits
include("WalkAndInterp.jl")
include("ComputePressurisedCrackDn.jl")
include("anneal.jl")
include("PropagateFracture.jl")
include("StrainEnergyRelease.jl")
include("FindPropAngleAndPoint.jl")
println("Scale length in here relative to energy")

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
include("PollardSegall1987_FractureSlipProfile.jl")
include("Burgmann1994_FractureLinearFrictionSlipProfile.jl")


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