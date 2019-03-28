#Setup the traction vector boundary condition depending on the input type 
function SetupTractionVector(BoundaryConditions::Stresses,FaceNormalVector,λ,G,n)

	if typeof(BoundaryConditions.σxx) == Float64 || Int64
		RepeatStruct(BoundaryConditions,n)
	end

	( Tn,Tds,Tss ) = CalculateNormalAndShearTractions3d( BoundaryConditions,FaceNormalVector);
	return Tn, Tds, Tss

end

function SetupTractionVector(BoundaryConditions::Strains,FaceNormalVector,λ,G,n)

	if typeof(BoundaryConditions.εxx) == Float64 || Int64
		RepeatStruct(BoundaryConditions,n)
	end

	#Convert to stress
	(StressTensor) = CutAndDisplaceJulia.HookesLaw3dStrain2Stress(BoundaryConditions,FaceNormalVector,λ,G);
	( Tn,Tds,Tss )=SetupTractionVector(BoundaryConditions)
	return Tn, Tds, Tss

end

function SetupTractionVector(BoundaryConditions::Tractions,FaceNormalVector,λ,G,n)

	if typeof(BoundaryConditions.Tn) == Float64 || Int64
		RepeatStruct(BoundaryConditions,n)
	end

	Tn=BoundaryConditions.Tn;
	Tss=BoundaryConditions.Tss;
	Tds=BoundaryConditions.Tds;
	return Tn, Tds, Tss

end

function SetupTractionVector(BoundaryConditions::MixedBoundaryConditions,FaceNormalVector,λ,G,n)

	(σTn,σTds,σTss)=SetupTractionVector(BoundaryConditions.Stresses)
	(Tn,Tds,Tss)=SetupTractionVector(BoundaryConditions.Tractions)
	Tn=	Tn.+σTn;
	Tds=Tds.+σTds;	
	Tss=Tss.+σTss;

	return Tn, Tds, Tss

end