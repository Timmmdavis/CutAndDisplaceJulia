#If we are not looking at fixed els
function CalculateInfluenceMatrices3D(FaceNormalVector,MidPoint,P1,P2,P3,ν,G,λ,FixedEls,HSFlag,n)

	(x,y,z,DssVec,DdsVec,DnVec,StressFlag,CosAx,CosAy,CosAz)=SetupCollationPoints(FaceNormalVector,MidPoint,n)

	println("Calling TD")
	#Computing for fixed els
	FixedFlag=FixedEls.==1
	NotFixedFlag=FixedEls.!=1;
	if any(FixedFlag)

		#Compute displacement on fixed elements
		xFix=x[FixedFlag];
		yFix=y[FixedFlag];
		zFix=z[FixedFlag];
		DispFlag=1;
		StrainFlag=0;
		(StrainInfMat,DispInfMat)=TD(xFix,yFix,zFix,P1,P2,P3,DssVec,DdsVec,DnVec,ν,G,DispFlag,StrainFlag,HSFlag)		
		#Compute stress on non fixed elements
		xNoFix=x[NotFixedFlag];
		yNoFix=y[NotFixedFlag];
		zNoFix=z[NotFixedFlag];	
		DispFlag=0;

		StrainFlag=1;
		(StrainInfMat,~)=TD(xNoFix,yNoFix,zNoFix,P1,P2,P3,DssVec,DdsVec,DnVec,ν,G,DispFlag,StrainFlag,HSFlag);
		 
		TractionInfMats=ConvertStrainInfMatsToTraction(StrainInfMat,λ,G,CosAx,CosAy,CosAz);

		#Now Putting DispInfMats into TractionInfMats
		TractionInfMats=ConcatInfMats(TractionInfMats,DispInfMat,NotFixedFlag,FixedFlag,n)

		return TractionInfMats 

	else

		DispFlag=0;
		#Call TDE function
		(StrainInfMat,DispInfMat)=TD(x,y,z,P1,P2,P3,DssVec,DdsVec,DnVec,ν,G,DispFlag,StressFlag,HSFlag)

		 TractionInfMats=ConvertStrainInfMatsToTraction(StrainInfMat,λ,G,CosAx,CosAy,CosAz)

		return TractionInfMats 

	end


end



function SetupCollationPoints(FaceNormalVector,MidPoint,n)
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
	#Flags
	StressFlag=1;
	return x,y,z,DssVec,DdsVec,DnVec,StressFlag,CosAx,CosAy,CosAz
end

function ConvertStrainInfMatsToTraction(SIM,λ,G,CosAx,CosAy,CosAz)

	#Faster threaded version for large mats. 

	DssTn=zeros(size(SIM.εxxDss));DdsTn=copy(DssTn);	DnTn=copy(DssTn);
	DssTss=copy(DssTn);		DdsTss=copy(DssTn);		DnTss=copy(DssTn);
	DssTds=copy(DssTn);		DdsTds=copy(DssTn);		DnTds=copy(DssTn);
	
	#Calculates the directions of the dipslip and ss directions
	(StrikeSlipCosine,DipSlipCosine) = CalculateSSandDSDirs( CosAx,CosAy,CosAz );

	Threads.@threads for i=1:size(SIM.εxxDss,2) 

	#Safe for threading 
	Tx=zeros(size(SIM.εxxDss,1));
	Ty=zeros(size(SIM.εxxDss,1));
	Tz=zeros(size(SIM.εxxDss,1));

	#Converting strains to stress tensor influences  
	(σxxDss,σyyDss,σzzDss,σxyDss,σxzDss,σyzDss) = HookesLaw3DStrain2Stress!(view(SIM.εxxDss,:,i),view(SIM.εyyDss,:,i),view(SIM.εzzDss,:,i),view(SIM.εxyDss,:,i),view(SIM.εxzDss,:,i),view(SIM.εyzDss,:,i),λ,G);
	(σxxDds,σyyDds,σzzDds,σxyDds,σxzDds,σyzDds) = HookesLaw3DStrain2Stress!(view(SIM.εxxDds,:,i),view(SIM.εyyDds,:,i),view(SIM.εzzDds,:,i),view(SIM.εxyDds,:,i),view(SIM.εxzDds,:,i),view(SIM.εyzDds,:,i),λ,G);
	(σxxDn,σyyDn,σzzDn,σxyDn,σxzDn,σyzDn) 		= HookesLaw3DStrain2Stress!(view(SIM.εxxDn,:,i) ,view(SIM.εyyDn,:,i), view(SIM.εzzDn,:,i) ,view(SIM.εxyDn,:,i) ,view(SIM.εxzDn,:,i) ,view(SIM.εyzDn,:,i), λ,G);
	#Compute normal traction
	CalculateNormalTraction3D!( σxxDss,σyyDss,σzzDss,σxyDss,σxzDss,σyzDss,CosAx,CosAy,CosAz,view(DssTn,:,i))
	CalculateNormalTraction3D!( σxxDds,σyyDds,σzzDds,σxyDds,σxzDds,σyzDds,CosAx,CosAy,CosAz,view(DdsTn,:,i))
	CalculateNormalTraction3D!( σxxDn,σyyDn,σzzDn,σxyDn,σxzDn,σyzDn,CosAx,CosAy,CosAz,      view(DnTn,:,i) )

	TractionVectorCartesianComponents3D!(σxxDss,σyyDss,σzzDss,σxyDss,σxzDss,σyzDss,CosAx,CosAy,CosAz,Tx,Ty,Tz)
	CalculateTractionInChosenDirection3D!( Tx,Ty,Tz,CosAx,CosAy,CosAz,StrikeSlipCosine,view(DssTss,:,i) );
	CalculateTractionInChosenDirection3D!( Tx,Ty,Tz,CosAx,CosAy,CosAz,DipSlipCosine,   view(DssTds,:,i) );

	TractionVectorCartesianComponents3D!(σxxDds,σyyDds,σzzDds,σxyDds,σxzDds,σyzDds,CosAx,CosAy,CosAz,Tx,Ty,Tz)
	CalculateTractionInChosenDirection3D!( Tx,Ty,Tz,CosAx,CosAy,CosAz,StrikeSlipCosine,view(DdsTss,:,i) );
	CalculateTractionInChosenDirection3D!( Tx,Ty,Tz,CosAx,CosAy,CosAz,DipSlipCosine,   view(DdsTds,:,i) );

	TractionVectorCartesianComponents3D!(σxxDn,σyyDn,σzzDn,σxyDn,σxzDn,σyzDn,CosAx,CosAy,CosAz,Tx,Ty,Tz)
	CalculateTractionInChosenDirection3D!( Tx,Ty,Tz,CosAx,CosAy,CosAz,StrikeSlipCosine,view(DnTss,:,i) );
	CalculateTractionInChosenDirection3D!( Tx,Ty,Tz,CosAx,CosAy,CosAz,DipSlipCosine,   view(DnTds,:,i) );

	end
	#Now putting influence matricies inside a predefined structure
	TractionInfMats = TractionInf(	DnTn,DnTss,DnTds, 	DssTn,DssTss,DssTds, 	DdsTn,DdsTss,DdsTds);

#=
	#Converting strains to stress tensor influences  
	(σxxDss,σyyDss,σzzDss,σxyDss,σxzDss,σyzDss) = HookesLaw3DStrain2Stress(SIM.εxxDss,SIM.εyyDss,SIM.εzzDss,SIM.εxyDss,SIM.εxzDss,SIM.εyzDss,λ,G);
	(σxxDds,σyyDds,σzzDds,σxyDds,σxzDds,σyzDds) = HookesLaw3DStrain2Stress(SIM.εxxDds,SIM.εyyDds,SIM.εzzDds,SIM.εxyDds,SIM.εxzDds,SIM.εyzDds,λ,G);
	(σxxDn,σyyDn,σzzDn,σxyDn,σxzDn,σyzDn) 		= HookesLaw3DStrain2Stress(SIM.εxxDn,SIM.εyyDn,SIM.εzzDn,SIM.εxyDn,SIM.εxzDn,SIM.εyzDn,λ,G);
	#Compute normal traction
	DssTn=CalculateNormalTraction3D( σxxDss,σyyDss,σzzDss,σxyDss,σxzDss,σyzDss,CosAx,CosAy,CosAz )
	DdsTn=CalculateNormalTraction3D( σxxDds,σyyDds,σzzDds,σxyDds,σxzDds,σyzDds,CosAx,CosAy,CosAz )
	DnTn =CalculateNormalTraction3D( σxxDn,σyyDn,σzzDn,σxyDn,σxzDn,σyzDn,CosAx,CosAy,CosAz )
	#Cart components of traction vector
	(DssT1x,DssT2y,DssT3z ) =TractionVectorCartesianComponents3D(σxxDss,σyyDss,σzzDss,σxyDss,σxzDss,σyzDss,CosAx,CosAy,CosAz)
	(DdsT1x,DdsT2y,DdsT3z ) =TractionVectorCartesianComponents3D(σxxDds,σyyDds,σzzDds,σxyDds,σxzDds,σyzDds,CosAx,CosAy,CosAz)
	(DnT1x,DnT2y,DnT3z ) 	=TractionVectorCartesianComponents3D(σxxDn,σyyDn,σzzDn,σxyDn,σxzDn,σyzDn,CosAx,CosAy,CosAz)
	#Calculates the directions of the dipslip and ss directions
	(StrikeSlipCosine,DipSlipCosine) = CalculateSSandDSDirs( CosAx,CosAy,CosAz );
	#Compute strike slip traction
	( DssTss ) = CalculateTractionInChosenDirection3D( DssT1x,DssT2y,DssT3z,CosAx,CosAy,CosAz,StrikeSlipCosine );
	( DdsTss ) = CalculateTractionInChosenDirection3D( DdsT1x,DdsT2y,DdsT3z,CosAx,CosAy,CosAz,StrikeSlipCosine );
	( DnTss  ) = CalculateTractionInChosenDirection3D( DnT1x,DnT2y,DnT3z,CosAx,CosAy,CosAz,StrikeSlipCosine );
	#Compute dip slip traction
	( DssTds ) = CalculateTractionInChosenDirection3D( DssT1x,DssT2y,DssT3z,CosAx,CosAy,CosAz,DipSlipCosine );
	( DdsTds ) = CalculateTractionInChosenDirection3D( DdsT1x,DdsT2y,DdsT3z,CosAx,CosAy,CosAz,DipSlipCosine );
	( DnTds  ) = CalculateTractionInChosenDirection3D( DnT1x,DnT2y,DnT3z,CosAx,CosAy,CosAz,DipSlipCosine );
	#Now putting influence matricies inside a predefined structure
	TractionInfMats = TractionInf( DnTn,DnTss,DnTds, DssTn,DssTss,DssTds, DdsTn,DdsTss,DdsTds);
=#
end


function PutDispInfsIntoInfMat(TractionInfMat,DispInfMat,TractionFlag,DispFlag,n)
#Fills disp rows in correct place with the traction influence matrix.
	
	#This is allocated here. Pre allocated INF mats could be passed into TD if we were concerned with speed here. 
	FullInfMat=zeros(n,n);

	count=0;
	for j=1:length(DispFlag)
		if DispFlag[j]==true	
			count+=1;
			for i=1:n
				FullInfMat[j,i]=DispInfMat[count,i]
			end
		end
	end
	count=0;
	for j=1:length(TractionFlag)
		if TractionFlag[j]==true		
			count+=1;		
			for i=1:n
				FullInfMat[j,i]=TractionInfMat[count,i]
			end
		end
	end	
	#Needs to be a single vector
	TractionFlag=reshape(TractionFlag,length(TractionFlag))
	#Now remove cols that we dont need
	FullInfMat=FullInfMat[:,TractionFlag];

	return FullInfMat

end


#Repeat values in structures to match a given size (assuming your structures are mutable)
#This assumes you only want to repeat the array along a single dimension
function ConcatInfMats(TractionStruc,DispStruct,TractionFlag,DispFlag,n)

	#Get fields in structure
	TractionFields=fieldnames(typeof(TractionStruc));
	DispFields=fieldnames(typeof(DispStruct));

	
	for i=1:length(TractionFields)
		#Check field i
		TractionMat=getfield(TractionStruc, TractionFields[i])
		DispMat=getfield(DispStruct, DispFields[i])

		TractionMat=PutDispInfsIntoInfMat(TractionMat,DispMat,TractionFlag,DispFlag,n)
		setfield!(TractionStruc,TractionFields[i],TractionMat)
	end

	return TractionStruc

end