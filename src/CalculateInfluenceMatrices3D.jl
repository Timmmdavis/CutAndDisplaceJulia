#If we are not looking at fixed els
function CalculateInfluenceMatrices3D(FaceNormalVector,MidPoint,P1,P2,P3,ν,G,λ,FixedEls,HSFlag,n)

	(x,y,z,DssVec,DdsVec,DnVec,StressFlag,CosAx,CosAy,CosAz)=SetupCollationPoints(FaceNormalVector,MidPoint,n)
	
	#Computing for fixed els
	FixedFlag=FixedEls.==1
	NotFixedFlag=FixedEls.!=1;
	if any(FixedFlag)


		#Compute displacement on fixed elements
		xFix=x[FixedFlag];
		yFix=y[FixedFlag];
		zFix=z[FixedFlag];
		xFix=reshape(xFix,length(xFix),1)
		yFix=reshape(yFix,length(yFix),1)
		zFix=reshape(zFix,length(zFix),1)

		DispFlag=1;
		StressFlag=0;
		(εxxDn,εyyDn,εzzDn,εxyDn,εxzDn,εyzDn,
		 εxxDss,εyyDss,εzzDss,εxyDss,εxzDss,εyzDss,
		 εxxDds,εyyDds,εzzDds,εxyDds,εxzDds,εyzDds,
		 DnUx,DnUy,DnUz,
		 DssUx,DssUy,DssUz,
		 DdsUx,DdsUy,DdsUz)=TD(xFix,yFix,zFix,P1,P2,P3,DssVec,DdsVec,DnVec,ν,G,DispFlag,StressFlag,HSFlag)
		#Compute stress on non fixed elements
		xNoFix=x[NotFixedFlag];
		yNoFix=y[NotFixedFlag];
		zNoFix=z[NotFixedFlag];	
		xNoFix=reshape(xNoFix,length(xNoFix),1)
		yNoFix=reshape(yNoFix,length(yNoFix),1)
		zNoFix=reshape(zNoFix,length(zNoFix),1)

		DispFlag=0;
		StressFlag=1;
		(εxxDn,εyyDn,εzzDn,εxyDn,εxzDn,εyzDn,
		 εxxDss,εyyDss,εzzDss,εxyDss,εxzDss,εyzDss,
		 εxxDds,εyyDds,εzzDds,εxyDds,εxzDds,εyzDds)=
		 TD(xNoFix,yNoFix,zNoFix,P1,P2,P3,DssVec,DdsVec,DnVec,ν,G,DispFlag,StressFlag,HSFlag);
		 
		TractionInfMats=ConvertInfMatsToTraction(εxxDn,εyyDn,εzzDn,εxyDn,εxzDn,εyzDn,
				 								  εxxDss,εyyDss,εzzDss,εxyDss,εxzDss,εyzDss,
				 								  εxxDds,εyyDds,εzzDds,εxyDds,εxzDds,εyzDds,
				 								  λ,G,CosAx,CosAy,CosAz);

		DispInfMats = DispInf(DnUx,DnUy,DnUz,DssUx,DssUy,DssUz,DdsUx,DdsUy,DdsUz);

		#Now Putting DispInfMats into TractionInfMats
		TractionInfMats=ConcatInfMats(TractionInfMats,DispInfMats,NotFixedFlag,FixedFlag,n)

		return TractionInfMats 

	else

		DispFlag=0;
		#Call TDE function
		(εxxDn,εyyDn,εzzDn,εxyDn,εxzDn,εyzDn,
		 εxxDss,εyyDss,εzzDss,εxyDss,εxzDss,εyzDss,
		 εxxDds,εyyDds,εzzDds,εxyDds,εxzDds,εyzDds)= TD(x,y,z,P1,P2,P3,DssVec,DdsVec,DnVec,ν,G,DispFlag,StressFlag,HSFlag)
		 
		 TractionInfMats=ConvertInfMatsToTraction(εxxDn,εyyDn,εzzDn,εxyDn,εxzDn,εyzDn,
				 								  εxxDss,εyyDss,εzzDss,εxyDss,εxzDss,εyzDss,
				 								  εxxDds,εyyDds,εzzDds,εxyDds,εxzDds,εyzDds,
				 								  λ,G,CosAx,CosAy,CosAz)

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

function ConvertInfMatsToTraction(εxxDn,εyyDn,εzzDn,εxyDn,εxzDn,εyzDn,
								  εxxDss,εyyDss,εzzDss,εxyDss,εxzDss,εyzDss,
								  εxxDds,εyyDds,εzzDds,εxyDds,εxzDds,εyzDds,
								  λ,G,CosAx,CosAy,CosAz)
	#Converting strains to stress tensor influences  
	(σxxDss,σyyDss,σzzDss,σxyDss,σxzDss,σyzDss) = HookesLaw3dStrain2Stress(εxxDss,εyyDss,εzzDss,εxyDss,εxzDss,εyzDss,λ,G);
	(σxxDds,σyyDds,σzzDds,σxyDds,σxzDds,σyzDds) = HookesLaw3dStrain2Stress(εxxDds,εyyDds,εzzDds,εxyDds,εxzDds,εyzDds,λ,G);
	(σxxDn,σyyDn,σzzDn,σxyDn,σxzDn,σyzDn) 		= HookesLaw3dStrain2Stress(εxxDn,εyyDn,εzzDn,εxyDn,εxzDn,εyzDn,λ,G);
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
	TractionInfMats = TractionInf(	DnTn,DnTss,DnTds, 	DssTn,DssTss,DssTds, 	DdsTn,DdsTss,DdsTds);
end


function PutDispInfsIntoInfMat(TractionInfMat,DispInfMat,TractionFlag,DispFlag,n)
#Fills disp rows in correct place.

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