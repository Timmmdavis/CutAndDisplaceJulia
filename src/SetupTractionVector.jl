#Setup the traction vector boundary condition depending on the input type 
function SetupTractionVector(BC::Stresses,FaceNormalVector,n,G,λ)
							#Importing BoundaryConditions as BC

	RepeatStruct(BC,n)
	( Tn,Tds,Tss ) = CalculateNormalAndShearTractions3D( BC.σxx,BC.σyy,BC.σzz,BC.σxy,BC.σxz,BC.σyz,FaceNormalVector);
	return Tn, Tds, Tss

end

function SetupTractionVector(BoundaryConditions::Strains,FaceNormalVector,n,G,λ)

	RepeatStruct(BoundaryConditions,n)

	#Convert to stress
	(StressTensor) = HookesLaw3DStrain2Stress(BoundaryConditions,FaceNormalVector,λ,G);
	( Tn,Tds,Tss )=SetupTractionVector(BoundaryConditions)
	return Tn, Tds, Tss

end

function SetupTractionVector(BoundaryConditions::Tractions,FaceNormalVector,n,G,λ)

	RepeatStruct(BoundaryConditions,n)

	Tn=BoundaryConditions.Tn;
	Tss=BoundaryConditions.Tss;
	Tds=BoundaryConditions.Tds;
	return Tn, Tds, Tss

end

function SetupTractionVector(BoundaryConditions::MixedBoundaryConditions,FaceNormalVector,n,G,λ)

	(σTn,σTds,σTss)=SetupTractionVector(BoundaryConditions.Stresses)
	(Tn,Tds,Tss)=SetupTractionVector(BoundaryConditions.Tractions)
	Tn=	Tn.+σTn;
	Tds=Tds.+σTds;	
	Tss=Tss.+σTss;

	return Tn, Tds, Tss

end