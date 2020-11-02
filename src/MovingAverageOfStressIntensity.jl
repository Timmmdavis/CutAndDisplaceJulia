function MovingAverageOfStressIntensity(avgeverynth,P1,P2,P3,FaceNormalVector,MidPoint,FeP1P2S,FeP1P3S,FeP2P3S)

#avgeverynth - if 3 averages from the surrounding triangles
#		 if 5 averages from two tris either side


if avgeverynth==3
	(UniqueEdges,LeadingPoints,TrailingPoints,InnerPoints,rerunFunc,P1,P2,P3,FaceNormalVector,MidPoint)=
	CutAndDisplaceJulia.GetSortedEdgesOfMeshList(P1,P2,P3,FaceNormalVector,MidPoint)
	(SortedTriangles,ConnectedEdge)=CutAndDisplaceJulia.ConnectedConstraints(P1,P2,P3,MidPoint);
	LeadingPoint=[0. 0. 0.]
	TrailingPoint=[0. 0. 0.]
	InnerPoint   =[0. 0. 0.]    
	LeadingPointOld  =[NaN NaN NaN] 
	TrailingPointOld =[NaN NaN NaN] 
	InnerPointOld    =[NaN NaN NaN] 
	BackPoint    =[NaN NaN NaN] 
	lps=0

	FeP2P3S_StrainEnergy=copy(FeP2P3S.StrainEnergy)
	FeP2P3S_K1=copy(FeP2P3S.K1)
	FeP2P3S_K2=copy(FeP2P3S.K2)
	FeP2P3S_K3=copy(FeP2P3S.K3)

	FeP1P3S_StrainEnergy=copy(FeP1P3S.StrainEnergy)
	FeP1P3S_K1=copy(FeP1P3S.K1)
	FeP1P3S_K2=copy(FeP1P3S.K2)
	FeP1P3S_K3=copy(FeP1P3S.K3)

	FeP1P2S_StrainEnergy=copy(FeP1P2S.StrainEnergy)
	FeP1P2S_K1=copy(FeP1P2S.K1)	
	FeP1P2S_K2=copy(FeP1P2S.K2)	
	FeP1P2S_K3=copy(FeP1P2S.K3)	
	#For each edge loop
	for i=1:length(UniqueEdges)
	  b=vec(UniqueEdges[i])
	  for j=b
	  	#Get the two triangles surrounding j
	  	j_trailing=j-1
	  	j_current=j
	  	j_future=j+1
	  	if j==1
	  		j_trailing=b[end]
	  	elseif j==b[end]
	  		j_future=b[1]
	  	end
	    #Extract the points on the current bit of the edge
	    (k,Idx_trailing) =GrabPointNew6(InnerPoints,P1,P2,P3,j_trailing)
	    if k==1
	    	Vtrailing=FeP2P3S_StrainEnergy[Idx_trailing]
	    	VtrailingK1=FeP2P3S_K1[Idx_trailing]
	    	VtrailingK2=FeP2P3S_K2[Idx_trailing]
	    	VtrailingK3=FeP2P3S_K3[Idx_trailing]
	    elseif k==2
	    	Vtrailing=FeP1P3S_StrainEnergy[Idx_trailing]
	    	VtrailingK1=FeP1P3S_K1[Idx_traili_g]
	    	VtrailingK2=FeP1P3S_K2[Idx_trailing]
	    	VtrailingK3=FeP1P3S_K3[Idx_trailing]	    	
	    elseif k==3
	    	Vtrailing=FeP1P2S_StrainEnergy[Idx_trailing]
	    	VtrailingK1=FeP1P2S_K1[Idx_trailing]
	    	VtrailingK2=FeP1P2S_K2[Idx_trailing]
	    	VtrailingK3=FeP1P2S_K3[Idx_trailing]		    	
	    end    	
	    (k,Idx_future) =GrabPointNew6(InnerPoints,P1,P2,P3,j_future)
		if k==1
	    	Vfuture=FeP2P3S_StrainEnergy[Idx_future]
	    	VfutureK1=FeP2P3S_K1[Idx_future]
	    	VfutureK2=FeP2P3S_K2[Idx_future]
	    	VfutureK3=FeP2P3S_K3[Idx_future]	    	
	    elseif k==2
	    	Vfuture=FeP1P3S_StrainEnergy[Idx_future]
	    	VfutureK1=FeP1P3S_K1[Idx_future]
	    	VfutureK2=FeP1P3S_K2[Idx_future]
	    	VfutureK3=FeP1P3S_K3[Idx_future]		    	
	    elseif k==3
	    	Vfuture=FeP1P2S_StrainEnergy[Idx_future]
	    	VfutureK1=FeP1P2S_K1[Idx_future]
	    	VfutureK2=FeP1P2S_K2[Idx_future]
	    	VfutureK3=FeP1P2S_K3[Idx_future]		    	
	    end   
	    (k,Idx_current) =GrabPointNew6(InnerPoints,P1,P2,P3,j_current)
		if k==1
			avg=(Vtrailing+FeP2P3S_StrainEnergy[Idx_current]+Vfuture)/3;
			avgK1=(VtrailingK1+FeP2P3S_K1[Idx_current]+VfutureK1)/3;
			avgK2=(VtrailingK2+FeP2P3S_K2[Idx_current]+VfutureK2)/3;
			avgK3=(VtrailingK3+FeP2P3S_K3[Idx_current]+VfutureK3)/3;
			if isnan(avg)
				continue
			end
			FeP2P3S.StrainEnergy[Idx_current]=avg
			FeP2P3S.K1[Idx_current]=avgK1
			FeP2P3S.K2[Idx_current]=avgK2
			FeP2P3S.K3[Idx_current]=avgK3
	    elseif k==2
	    	avg=(Vtrailing+FeP1P3S_StrainEnergy[Idx_current]+Vfuture)/3
			avgK1=(VtrailingK1+FeP1P3S_K1[Idx_current]+VfutureK1)/3;
			avgK2=(VtrailingK2+FeP1P3S_K2[Idx_current]+VfutureK2)/3;
			avgK3=(VtrailingK3+FeP1P3S_K3[Idx_current]+VfutureK3)/3;	    	
			if isnan(avg)
				continue
			end
	    	FeP1P3S.StrainEnergy[Idx_current]=avg
			FeP1P3S.K1[Idx_current]=avgK1
			FeP1P3S.K2[Idx_current]=avgK2
			FeP1P3S.K3[Idx_current]=avgK3	    	
	    elseif k==3
	    	avg=(Vtrailing+FeP1P2S_StrainEnergy[Idx_current]+Vfuture)/3 
	    	avgK1=(VtrailingK1+FeP1P2S_K1[Idx_current]+VfutureK1)/3;
			avgK2=(VtrailingK2+FeP1P2S_K2[Idx_current]+VfutureK2)/3;
			avgK3=(VtrailingK3+FeP1P2S_K3[Idx_current]+VfutureK3)/3;
	    	if isnan(avg)
				continue
			end
	    	FeP1P2S.StrainEnergy[Idx_current]=avg  
	    	FeP1P2S.K1[Idx_current]=avgK1
			FeP1P2S.K2[Idx_current]=avgK2
			FeP1P2S.K3[Idx_current]=avgK3	
	    end       
	    #@info Idx_trailing Idx_current Idx_future
	    end
	end 
end

if avgeverynth==5
	(UniqueEdges,LeadingPoints,TrailingPoints,InnerPoints,rerunFunc,P1,P2,P3,FaceNormalVector,MidPoint)=
	CutAndDisplaceJulia.GetSortedEdgesOfMeshList(P1,P2,P3,FaceNormalVector,MidPoint)
	(SortedTriangles,ConnectedEdge)=CutAndDisplaceJulia.ConnectedConstraints(P1,P2,P3,MidPoint);
	LeadingPoint=[0. 0. 0.]
	TrailingPoint=[0. 0. 0.]
	InnerPoint   =[0. 0. 0.]    
	LeadingPointOld  =[NaN NaN NaN] 
	TrailingPointOld =[NaN NaN NaN] 
	InnerPointOld    =[NaN NaN NaN] 
	BackPoint    =[NaN NaN NaN] 
	lps=0

	FeP2P3S_StrainEnergy=copy(FeP2P3S.StrainEnergy)
	FeP2P3S_K1=copy(FeP2P3S.K1)
	FeP2P3S_K2=copy(FeP2P3S.K2)
	FeP2P3S_K3=copy(FeP2P3S.K3)

	FeP1P3S_StrainEnergy=copy(FeP1P3S.StrainEnergy)
	FeP1P3S_K1=copy(FeP1P3S.K1)
	FeP1P3S_K2=copy(FeP1P3S.K2)
	FeP1P3S_K3=copy(FeP1P3S.K3)

	FeP1P2S_StrainEnergy=copy(FeP1P2S.StrainEnergy)
	FeP1P2S_K1=copy(FeP1P2S.K1)	
	FeP1P2S_K2=copy(FeP1P2S.K2)	
	FeP1P2S_K3=copy(FeP1P2S.K3)	

	#For each edge loop
	for i=1:length(UniqueEdges)
	  b=vec(UniqueEdges[i])
	  for j=b
	  	#Get the two triangles surrounding j
	  	j_trailing2=j-2
	  	j_trailing=j-1
	  	j_current=j
	  	j_future=j+1
	  	j_future2=j+2
	  	if j==1
	  		j_trailing=b[end]
	  		j_trailing2=b[end-1]
	  	elseif j==b[end]
	  		j_future=b[1]
	  		j_future2=b[2]
	  	end
	  	if j==2
	  		j_trailing2=b[end]
	  	elseif j==b[end-1]
	  		j_future2=b[1]
	  	end

	    #Extract the points on the current bit of the edge
	    (k,Idx_trailing) =GrabPointNew6(InnerPoints,P1,P2,P3,j_trailing)
	    if k==1
	    	Vtrailing=FeP2P3S_StrainEnergy[Idx_trailing]
	    	VtrailingK1=FeP2P3S_K1[Idx_trailing]
	    	VtrailingK2=FeP2P3S_K2[Idx_trailing]
	    	VtrailingK3=FeP2P3S_K3[Idx_trailing]
	    elseif k==2
	    	Vtrailing=FeP1P3S_StrainEnergy[Idx_trailing]
	    	VtrailingK1=FeP1P3S_K1[Idx_traili_g]
	    	VtrailingK2=FeP1P3S_K2[Idx_trailing]
	    	VtrailingK3=FeP1P3S_K3[Idx_trailing]	    	
	    elseif k==3
	    	Vtrailing=FeP1P2S_StrainEnergy[Idx_trailing]
	    	VtrailingK1=FeP1P2S_K1[Idx_trailing]
	    	VtrailingK2=FeP1P2S_K2[Idx_trailing]
	    	VtrailingK3=FeP1P2S_K3[Idx_trailing]
	    end   
		(k,Idx_trailing) =GrabPointNew6(InnerPoints,P1,P2,P3,j_trailing2)
	    if k==1
	    	Vtrailing_2=FeP2P3S_StrainEnergy[Idx_trailing]
	    	VtrailingK1_2=FeP2P3S_K1[Idx_trailing]
	    	VtrailingK2_2=FeP2P3S_K2[Idx_trailing]
	    	VtrailingK3_2=FeP2P3S_K3[Idx_trailing]
	    elseif k==2
	    	Vtrailing_2=FeP1P3S_StrainEnergy[Idx_trailing]
	    	VtrailingK1_2=FeP1P3S_K1[Idx_traili_g]
	    	VtrailingK2_2=FeP1P3S_K2[Idx_trailing]
	    	VtrailingK3_2=FeP1P3S_K3[Idx_trailing]	    	
	    elseif k==3
	    	Vtrailing_2=FeP1P2S_StrainEnergy[Idx_trailing]
	    	VtrailingK1_2=FeP1P2S_K1[Idx_trailing]
	    	VtrailingK2_2=FeP1P2S_K2[Idx_trailing]
	    	VtrailingK3_2=FeP1P2S_K3[Idx_trailing]
	    end   
	    (k,Idx_future) =GrabPointNew6(InnerPoints,P1,P2,P3,j_future)
		if k==1
			Vfuture=FeP2P3S_StrainEnergy[Idx_future]
	    	VfutureK1=FeP2P3S_K1[Idx_future]
	    	VfutureK2=FeP2P3S_K2[Idx_future]
	    	VfutureK3=FeP2P3S_K3[Idx_future]	    	
	    elseif k==2
	    	Vfuture=FeP1P3S_StrainEnergy[Idx_future]
	    	VfutureK1=FeP1P3S_K1[Idx_future]
	    	VfutureK2=FeP1P3S_K2[Idx_future]
	    	VfutureK3=FeP1P3S_K3[Idx_future]		    	
	    elseif k==3
	    	Vfuture=FeP1P2S_StrainEnergy[Idx_future]
	    	VfutureK1=FeP1P2S_K1[Idx_future]
	    	VfutureK2=FeP1P2S_K2[Idx_future]
	    	VfutureK3=FeP1P2S_K3[Idx_future]		    	
	    end   
	    (k,Idx_future) =GrabPointNew6(InnerPoints,P1,P2,P3,j_future2)
		if k==1
			Vfuture_2=FeP2P3S_StrainEnergy[Idx_future]
	    	VfutureK1_2=FeP2P3S_K1[Idx_future]
	    	VfutureK2_2=FeP2P3S_K2[Idx_future]
	    	VfutureK3_2=FeP2P3S_K3[Idx_future]	    	
	    elseif k==2
	    	Vfuture_2=FeP1P3S_StrainEnergy[Idx_future]
	    	VfutureK1_2=FeP1P3S_K1[Idx_future]
	    	VfutureK2_2=FeP1P3S_K2[Idx_future]
	    	VfutureK3_2=FeP1P3S_K3[Idx_future]		    	
	    elseif k==3
	    	Vfuture_2=FeP1P2S_StrainEnergy[Idx_future]
	    	VfutureK1_2=FeP1P2S_K1[Idx_future]
	    	VfutureK2_2=FeP1P2S_K2[Idx_future]
	    	VfutureK3_2=FeP1P2S_K3[Idx_future]	
	    end   	    
	    (k,Idx_current) =GrabPointNew6(InnerPoints,P1,P2,P3,j_current)
		if k==1
			avg=(Vtrailing_2+Vtrailing+FeP2P3S_StrainEnergy[Idx_current]+Vfuture+Vfuture_2)/5;
			avgK1=(VtrailingK1_2+VtrailingK1+FeP2P3S_K1[Idx_current]+VfutureK1+VfutureK1_2)/5;
			avgK2=(VtrailingK2_2+VtrailingK2+FeP2P3S_K2[Idx_current]+VfutureK2+VfutureK2_2)/5;
			avgK3=(VtrailingK3_2+VtrailingK3+FeP2P3S_K3[Idx_current]+VfutureK3+VfutureK3_2)/5;
			if isnan(avg)
				continue
			end
			FeP2P3S.StrainEnergy[Idx_current]=avg
			FeP2P3S.K1[Idx_current]=avgK1
			FeP2P3S.K2[Idx_current]=avgK2
			FeP2P3S.K3[Idx_current]=avgK3
	    elseif k==2
			avg=(Vtrailing_2+Vtrailing+FeP1P3S_StrainEnergy[Idx_current]+Vfuture+Vfuture_2)/5;
			avgK1=(VtrailingK1_2+VtrailingK1+FeP1P3S_K1[Idx_current]+VfutureK1+VfutureK1_2)/5;
			avgK2=(VtrailingK2_2+VtrailingK2+FeP1P3S_K2[Idx_current]+VfutureK2+VfutureK2_2)/5;
			avgK3=(VtrailingK3_2+VtrailingK3+FeP1P3S_K3[Idx_current]+VfutureK3+VfutureK3_2)/5;
			if isnan(avg)
				continue
			end
			FeP1P3S.StrainEnergy[Idx_current]=avg
			FeP1P3S.K1[Idx_current]=avgK1
			FeP1P3S.K2[Idx_current]=avgK2
			FeP1P3S.K3[Idx_current]=avgK3
	    elseif k==3
			avg=(Vtrailing_2+Vtrailing+FeP1P2S_StrainEnergy[Idx_current]+Vfuture+Vfuture_2)/5;
			avgK1=(VtrailingK1_2+VtrailingK1+FeP1P2S_K1[Idx_current]+VfutureK1+VfutureK1_2)/5;
			avgK2=(VtrailingK2_2+VtrailingK2+FeP1P2S_K2[Idx_current]+VfutureK2+VfutureK2_2)/5;
			avgK3=(VtrailingK3_2+VtrailingK3+FeP1P2S_K3[Idx_current]+VfutureK3+VfutureK3_2)/5;
			if isnan(avg)
				continue
			end
			FeP1P2S.StrainEnergy[Idx_current]=avg
			FeP1P2S.K1[Idx_current]=avgK1
			FeP1P2S.K2[Idx_current]=avgK2
			FeP1P2S.K3[Idx_current]=avgK3
	    end       
	    #@info Idx_trailing Idx_current Idx_future
	    end
	end 
end

return FeP1P2S,FeP1P3S,FeP2P3S

end

function GrabPointNew6(PointsIdxList,P1,P2,P3,j)
#Extract the points on the current bit of the edge
InnerPointNo=0;
Indx=0
for k=1:3
    Idx=PointsIdxList[j,k]
    if Idx==0
        continue
    elseif k==1
        InnerPointNo=k
        Indx=Idx
    elseif k==2
        InnerPointNo=k
        Indx=Idx
    elseif k==3
        InnerPointNo=k
        Indx=Idx
    end
end

return InnerPointNo,Indx
end