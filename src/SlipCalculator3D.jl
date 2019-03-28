
function SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,DispFlag,StressFlag,HSFlag,BoundaryConditions)

	#Compute the tractions acting on the crack
	n=length(FaceNormalVector[:,1]);
	(Tn,Tds,Tss)=SetupTractionVector(BoundaryConditions,FaceNormalVector,λ,G,n);

	# Start some vectors (spaced points)
	CosAx =  FaceNormalVector[:,1];  
	CosAy =  FaceNormalVector[:,2];   
	CosAz =  FaceNormalVector[:,3];  
	x=zeros(n,1);
	y=zeros(n,1);
	z=zeros(n,1); 
	x[:]=MidPoint[:,1]-(CosAx*1e-12);
	y[:]=MidPoint[:,2]-(CosAy*1e-12);
	z[:]=MidPoint[:,3]-(CosAz*1e-12);

	#Repeating array if you want to test with more tris (not the solution vs okada, just speed)
	DssVec	=ones(n,1); #const 
	DdsVec	=ones(n,1); #const 
	DnVec	=ones(n,1); #const 

	println("Vars created -> to TD func")

	(εxxDn,εyyDn,εzzDn,εxyDn,εxzDn,εyzDn,
	 εxxDss,εyyDss,εzzDss,εxyDss,εxzDss,εyzDss,
	 εxxDds,εyyDds,εzzDds,εxyDds,εxzDds,εyzDds,
	 DnUx,DnUy,DnUz,
	 DssUx,DssUy,DssUz,
	 DdsUx,DdsUy,DdsUz)=
	 CutAndDisplaceJulia.TD(x,y,z,P1,P2,P3,DssVec,DdsVec,DnVec,ν,G,DispFlag,StressFlag,HSFlag)
	 
	println("Out of TD func") 

	#Converting strains to stress tensor influences  
	@time (σxxDss,σyyDss,σzzDss,σxyDss,σxzDss,σyzDss) = CutAndDisplaceJulia.HookesLaw3dStrain2Stress(εxxDss,εyyDss,εzzDss,εxyDss,εxzDss,εyzDss,λ,G);
	(σxxDds,σyyDds,σzzDds,σxyDds,σxzDds,σyzDds) = CutAndDisplaceJulia.HookesLaw3dStrain2Stress(εxxDds,εyyDds,εzzDds,εxyDds,εxzDds,εyzDds,λ,G);
	(σxxDn,σyyDn,σzzDn,σxyDn,σxzDn,σyzDn) 		= CutAndDisplaceJulia.HookesLaw3dStrain2Stress(εxxDn,εyyDn,εzzDn,εxyDn,εxzDn,εyzDn,λ,G);

	#Compute normal traction
	DssTn=CutAndDisplaceJulia.CalculateNormalTraction3d( σxxDss,σyyDss,σzzDss,σxyDss,σxzDss,σyzDss,CosAx,CosAy,CosAz )
	DdsTn=CutAndDisplaceJulia.CalculateNormalTraction3d( σxxDds,σyyDds,σzzDds,σxyDds,σxzDds,σyzDds,CosAx,CosAy,CosAz )
	DnTn =CutAndDisplaceJulia.CalculateNormalTraction3d( σxxDn,σyyDn,σzzDn,σxyDn,σxzDn,σyzDn,CosAx,CosAy,CosAz )

	#Cart components of traction vector
	(DssT1x,DssT2y,DssT3z ) =CutAndDisplaceJulia.TractionVectorCartesianComponents3d(σxxDss,σyyDss,σzzDss,σxyDss,σxzDss,σyzDss,CosAx,CosAy,CosAz)
	(DdsT1x,DdsT2y,DdsT3z ) =CutAndDisplaceJulia.TractionVectorCartesianComponents3d(σxxDds,σyyDds,σzzDds,σxyDds,σxzDds,σyzDds,CosAx,CosAy,CosAz)
	(DnT1x,DnT2y,DnT3z ) 	=CutAndDisplaceJulia.TractionVectorCartesianComponents3d(σxxDn,σyyDn,σzzDn,σxyDn,σxzDn,σyzDn,CosAx,CosAy,CosAz)

	#Calculates the directions of the dipslip and ss directions
	(StrikeSlipCosine,DipSlipCosine) = CutAndDisplaceJulia.CalculateSSandDSDirs( CosAx,CosAy,CosAz );

	#Compute strike slip traction
	( DssTss ) = CutAndDisplaceJulia.CalculateTractionInChosenDirection3d( DssT1x,DssT2y,DssT3z,CosAx,CosAy,CosAz,StrikeSlipCosine );
	( DdsTss ) = CutAndDisplaceJulia.CalculateTractionInChosenDirection3d( DdsT1x,DdsT2y,DdsT3z,CosAx,CosAy,CosAz,StrikeSlipCosine );
	( DnTss  ) = CutAndDisplaceJulia.CalculateTractionInChosenDirection3d( DnT1x,DnT2y,DnT3z,CosAx,CosAy,CosAz,StrikeSlipCosine );

	#Compute dip slip traction
	( DssTds ) = CutAndDisplaceJulia.CalculateTractionInChosenDirection3d( DssT1x,DssT2y,DssT3z,CosAx,CosAy,CosAz,DipSlipCosine );
	( DdsTds ) = CutAndDisplaceJulia.CalculateTractionInChosenDirection3d( DdsT1x,DdsT2y,DdsT3z,CosAx,CosAy,CosAz,DipSlipCosine );
	( DnTds  ) = CutAndDisplaceJulia.CalculateTractionInChosenDirection3d( DnT1x,DnT2y,DnT3z,CosAx,CosAy,CosAz,DipSlipCosine );

	#Now putting influence matricies inside a predefined structure
	TractionInfMats = TractionInf(	DnTn,DnTss,DnTds, 	DssTn,DssTss,DssTds, 	DdsTn,DdsTss,DdsTds);
	DispInfMats   = DispInf(	DnUx,DnUy,DnUz,		DssUx,DssUy,DssUz,		DdsUx,DdsUy,DdsUz);

	#Concatenate influence matrix ready for equation 
	Atn  = [-TractionInfMats.DnTn  -TractionInfMats.DssTn  -TractionInfMats.DdsTn ];     
	Atss = [-TractionInfMats.DnTss -TractionInfMats.DssTss -TractionInfMats.DdsTss];     
	Atds = [-TractionInfMats.DnTds -TractionInfMats.DssTds -TractionInfMats.DdsTds];     
	A= [Atn;Atss;Atds];  

	#Prep traction vector
	b=[Tn; Tss; Tds];

	#Do linear equation
	D=A\b;

	#Extract arrays
	Dn=D[1:n];
	Dss=D[n+1:2*n];
	Dds=D[n*2+1:3*n];

	return Dn, Dss, Dds

end

#Setup the traction vector depending on the input type 
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


#Repeat values in structures to match a given size (assuming your structures are mutable)
#This assumes you only want to repeat the array along a single dimension
function RepeatStruct(Structure,Size)

	#Get fields in structure
	FieldsInStruct=fieldnames(typeof(Structure));
	
	for i=1:length(FieldsInStruct)
		#Check field i
		Value=getfield(Structure, FieldsInStruct[i])
		#If field is a float then we repeat it
		if typeof(Value) == Float64 || Int64
			#repeat to predefined size
			Value=repeat([Value],Size[1])
			#Put back in structure
			setfield!(Structure,FieldsInStruct[i],Value)
		else
		#Call the function again for inner structure
			RepeatStruct(Value,Size)
		end
	end

	return Structure

end