function MovingAverageOfStressIntensity(avgeverynth,P1,P2,P3,FaceNormalVector,MidPoint,FeP1P2S,FeP1P3S,FeP2P3S;UseEdgeLength=true)

#avgeverynth - if 3 averages from the surrounding triangles
#		 if 5 averages from two tris either side
#UseEdgeLength - true: weight each edge by its length (FeLe), false: plain mean

UseEdgeLength=false

if avgeverynth==3
	(UniqueEdges,LeadingPoints,TrailingPoints,InnerPoints,rerunFunc,P1,P2,P3,FaceNormalVector,MidPoint)=
	CutAndDisplaceJulia.GetSortedEdgesOfMeshList(P1,P2,P3,FaceNormalVector,MidPoint)
	if rerunFunc==1
		error("Mesh needed cleaning to walk the boundary, indices no longer match the Fe arrays")
	end
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
	  n=length(b)
	  if n<3
	  	#Loop shorter than the window, neighbours would be counted twice
	  	continue
	  end
	  for p=1:n
	  	#Get the two triangles surrounding b[p]
	  	#p is the position inside this loop, mod1 wraps round this loop only
	  	j_trailing=b[mod1(p-1,n)]
	  	j_current=b[p]
	  	j_future=b[mod1(p+1,n)]
	    #Extract the points on the current bit of the edge
	    (k,Idx_trailing) =GrabPointNew6(InnerPoints,P1,P2,P3,j_trailing)
	    Ltrailing=EdgeLength(k,Idx_trailing,UseEdgeLength,FeP1P2S,FeP1P3S,FeP2P3S)
	    if k==1
	    	Vtrailing=FeP2P3S_StrainEnergy[Idx_trailing]
	    	VtrailingK1=FeP2P3S_K1[Idx_trailing]
	    	VtrailingK2=FeP2P3S_K2[Idx_trailing]
	    	VtrailingK3=FeP2P3S_K3[Idx_trailing]
	    elseif k==2
	    	Vtrailing=FeP1P3S_StrainEnergy[Idx_trailing]
	    	VtrailingK1=FeP1P3S_K1[Idx_trailing]
	    	VtrailingK2=FeP1P3S_K2[Idx_trailing]
	    	VtrailingK3=FeP1P3S_K3[Idx_trailing]	    	
	    elseif k==3
	    	Vtrailing=FeP1P2S_StrainEnergy[Idx_trailing]
	    	VtrailingK1=FeP1P2S_K1[Idx_trailing]
	    	VtrailingK2=FeP1P2S_K2[Idx_trailing]
	    	VtrailingK3=FeP1P2S_K3[Idx_trailing]		    	
	    end    	
	    (k,Idx_future) =GrabPointNew6(InnerPoints,P1,P2,P3,j_future)
	    Lfuture=EdgeLength(k,Idx_future,UseEdgeLength,FeP1P2S,FeP1P3S,FeP2P3S)
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
	    Lcurrent=EdgeLength(k,Idx_current,UseEdgeLength,FeP1P2S,FeP1P3S,FeP2P3S)
		if k==1
			avg=NanAvg((FeP2P3S_StrainEnergy[Idx_current],Vtrailing,Vfuture),(Lcurrent,Ltrailing,Lfuture));
			avgK1=NanAvg((FeP2P3S_K1[Idx_current],VtrailingK1,VfutureK1),(Lcurrent,Ltrailing,Lfuture));
			avgK2=NanAvg((FeP2P3S_K2[Idx_current],VtrailingK2,VfutureK2),(Lcurrent,Ltrailing,Lfuture));
			avgK3=NanAvg((FeP2P3S_K3[Idx_current],VtrailingK3,VfutureK3),(Lcurrent,Ltrailing,Lfuture));
			if isnan(avg)
				continue
			end
			FeP2P3S.StrainEnergy[Idx_current]=avg
			FeP2P3S.K1[Idx_current]=avgK1
			FeP2P3S.K2[Idx_current]=avgK2
			FeP2P3S.K3[Idx_current]=avgK3
	    elseif k==2
	    	avg=NanAvg((FeP1P3S_StrainEnergy[Idx_current],Vtrailing,Vfuture),(Lcurrent,Ltrailing,Lfuture));
			avgK1=NanAvg((FeP1P3S_K1[Idx_current],VtrailingK1,VfutureK1),(Lcurrent,Ltrailing,Lfuture));
			avgK2=NanAvg((FeP1P3S_K2[Idx_current],VtrailingK2,VfutureK2),(Lcurrent,Ltrailing,Lfuture));
			avgK3=NanAvg((FeP1P3S_K3[Idx_current],VtrailingK3,VfutureK3),(Lcurrent,Ltrailing,Lfuture));	    	
			if isnan(avg)
				continue
			end
	    	FeP1P3S.StrainEnergy[Idx_current]=avg
			FeP1P3S.K1[Idx_current]=avgK1
			FeP1P3S.K2[Idx_current]=avgK2
			FeP1P3S.K3[Idx_current]=avgK3	    	
	    elseif k==3
	    	avg=NanAvg((FeP1P2S_StrainEnergy[Idx_current],Vtrailing,Vfuture),(Lcurrent,Ltrailing,Lfuture));
	    	avgK1=NanAvg((FeP1P2S_K1[Idx_current],VtrailingK1,VfutureK1),(Lcurrent,Ltrailing,Lfuture));
			avgK2=NanAvg((FeP1P2S_K2[Idx_current],VtrailingK2,VfutureK2),(Lcurrent,Ltrailing,Lfuture));
			avgK3=NanAvg((FeP1P2S_K3[Idx_current],VtrailingK3,VfutureK3),(Lcurrent,Ltrailing,Lfuture));
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
	if rerunFunc==1
		error("Mesh needed cleaning to walk the boundary, indices no longer match the Fe arrays")
	end
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
	  n=length(b)
	  if n<5
	  	#Loop shorter than the window, neighbours would be counted twice
	  	continue
	  end
	  for p=1:n
	  	#Get the four triangles surrounding b[p]
	  	#p is the position inside this loop, mod1 wraps round this loop only
	  	j_trailing2=b[mod1(p-2,n)]
	  	j_trailing=b[mod1(p-1,n)]
	  	j_current=b[p]
	  	j_future=b[mod1(p+1,n)]
	  	j_future2=b[mod1(p+2,n)]

	    #Extract the points on the current bit of the edge
	    (k,Idx_trailing) =GrabPointNew6(InnerPoints,P1,P2,P3,j_trailing)
	    Ltrailing=EdgeLength(k,Idx_trailing,UseEdgeLength,FeP1P2S,FeP1P3S,FeP2P3S)
	    if k==1
	    	Vtrailing=FeP2P3S_StrainEnergy[Idx_trailing]
	    	VtrailingK1=FeP2P3S_K1[Idx_trailing]
	    	VtrailingK2=FeP2P3S_K2[Idx_trailing]
	    	VtrailingK3=FeP2P3S_K3[Idx_trailing]
	    elseif k==2
	    	Vtrailing=FeP1P3S_StrainEnergy[Idx_trailing]
	    	VtrailingK1=FeP1P3S_K1[Idx_trailing]
	    	VtrailingK2=FeP1P3S_K2[Idx_trailing]
	    	VtrailingK3=FeP1P3S_K3[Idx_trailing]	    	
	    elseif k==3
	    	Vtrailing=FeP1P2S_StrainEnergy[Idx_trailing]
	    	VtrailingK1=FeP1P2S_K1[Idx_trailing]
	    	VtrailingK2=FeP1P2S_K2[Idx_trailing]
	    	VtrailingK3=FeP1P2S_K3[Idx_trailing]
	    end   
		(k,Idx_trailing) =GrabPointNew6(InnerPoints,P1,P2,P3,j_trailing2)
		Ltrailing_2=EdgeLength(k,Idx_trailing,UseEdgeLength,FeP1P2S,FeP1P3S,FeP2P3S)
	    if k==1
	    	Vtrailing_2=FeP2P3S_StrainEnergy[Idx_trailing]
	    	VtrailingK1_2=FeP2P3S_K1[Idx_trailing]
	    	VtrailingK2_2=FeP2P3S_K2[Idx_trailing]
	    	VtrailingK3_2=FeP2P3S_K3[Idx_trailing]
	    elseif k==2
	    	Vtrailing_2=FeP1P3S_StrainEnergy[Idx_trailing]
	    	VtrailingK1_2=FeP1P3S_K1[Idx_trailing]
	    	VtrailingK2_2=FeP1P3S_K2[Idx_trailing]
	    	VtrailingK3_2=FeP1P3S_K3[Idx_trailing]	    	
	    elseif k==3
	    	Vtrailing_2=FeP1P2S_StrainEnergy[Idx_trailing]
	    	VtrailingK1_2=FeP1P2S_K1[Idx_trailing]
	    	VtrailingK2_2=FeP1P2S_K2[Idx_trailing]
	    	VtrailingK3_2=FeP1P2S_K3[Idx_trailing]
	    end   
	    (k,Idx_future) =GrabPointNew6(InnerPoints,P1,P2,P3,j_future)
	    Lfuture=EdgeLength(k,Idx_future,UseEdgeLength,FeP1P2S,FeP1P3S,FeP2P3S)
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
	    Lfuture_2=EdgeLength(k,Idx_future,UseEdgeLength,FeP1P2S,FeP1P3S,FeP2P3S)
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
	    #Don't reach across a gap: if the inner neighbour is NaN/Inf ignore the outer one too
	    if !isfinite(Vtrailing)
	    	Vtrailing_2=VtrailingK1_2=VtrailingK2_2=VtrailingK3_2=NaN
	    end
	    if !isfinite(Vfuture)
	    	Vfuture_2=VfutureK1_2=VfutureK2_2=VfutureK3_2=NaN
	    end
	    (k,Idx_current) =GrabPointNew6(InnerPoints,P1,P2,P3,j_current)
	    Lcurrent=EdgeLength(k,Idx_current,UseEdgeLength,FeP1P2S,FeP1P3S,FeP2P3S)
		if k==1
			avg=NanAvg((FeP2P3S_StrainEnergy[Idx_current],Vtrailing_2,Vtrailing,Vfuture,Vfuture_2),(Lcurrent,Ltrailing_2,Ltrailing,Lfuture,Lfuture_2));
			avgK1=NanAvg((FeP2P3S_K1[Idx_current],VtrailingK1_2,VtrailingK1,VfutureK1,VfutureK1_2),(Lcurrent,Ltrailing_2,Ltrailing,Lfuture,Lfuture_2));
			avgK2=NanAvg((FeP2P3S_K2[Idx_current],VtrailingK2_2,VtrailingK2,VfutureK2,VfutureK2_2),(Lcurrent,Ltrailing_2,Ltrailing,Lfuture,Lfuture_2));
			avgK3=NanAvg((FeP2P3S_K3[Idx_current],VtrailingK3_2,VtrailingK3,VfutureK3,VfutureK3_2),(Lcurrent,Ltrailing_2,Ltrailing,Lfuture,Lfuture_2));
			if isnan(avg)
				continue
			end
			FeP2P3S.StrainEnergy[Idx_current]=avg
			FeP2P3S.K1[Idx_current]=avgK1
			FeP2P3S.K2[Idx_current]=avgK2
			FeP2P3S.K3[Idx_current]=avgK3
	    elseif k==2
			avg=NanAvg((FeP1P3S_StrainEnergy[Idx_current],Vtrailing_2,Vtrailing,Vfuture,Vfuture_2),(Lcurrent,Ltrailing_2,Ltrailing,Lfuture,Lfuture_2));
			avgK1=NanAvg((FeP1P3S_K1[Idx_current],VtrailingK1_2,VtrailingK1,VfutureK1,VfutureK1_2),(Lcurrent,Ltrailing_2,Ltrailing,Lfuture,Lfuture_2));
			avgK2=NanAvg((FeP1P3S_K2[Idx_current],VtrailingK2_2,VtrailingK2,VfutureK2,VfutureK2_2),(Lcurrent,Ltrailing_2,Ltrailing,Lfuture,Lfuture_2));
			avgK3=NanAvg((FeP1P3S_K3[Idx_current],VtrailingK3_2,VtrailingK3,VfutureK3,VfutureK3_2),(Lcurrent,Ltrailing_2,Ltrailing,Lfuture,Lfuture_2));
			if isnan(avg)
				continue
			end
			FeP1P3S.StrainEnergy[Idx_current]=avg
			FeP1P3S.K1[Idx_current]=avgK1
			FeP1P3S.K2[Idx_current]=avgK2
			FeP1P3S.K3[Idx_current]=avgK3
	    elseif k==3
			avg=NanAvg((FeP1P2S_StrainEnergy[Idx_current],Vtrailing_2,Vtrailing,Vfuture,Vfuture_2),(Lcurrent,Ltrailing_2,Ltrailing,Lfuture,Lfuture_2));
			avgK1=NanAvg((FeP1P2S_K1[Idx_current],VtrailingK1_2,VtrailingK1,VfutureK1,VfutureK1_2),(Lcurrent,Ltrailing_2,Ltrailing,Lfuture,Lfuture_2));
			avgK2=NanAvg((FeP1P2S_K2[Idx_current],VtrailingK2_2,VtrailingK2,VfutureK2,VfutureK2_2),(Lcurrent,Ltrailing_2,Ltrailing,Lfuture,Lfuture_2));
			avgK3=NanAvg((FeP1P2S_K3[Idx_current],VtrailingK3_2,VtrailingK3,VfutureK3,VfutureK3_2),(Lcurrent,Ltrailing_2,Ltrailing,Lfuture,Lfuture_2));
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

function NanAvg(V,L)
#Length weighted mean of V[1] (the current edge) and whichever neighbours are not NaN.
#L holds the matching edge lengths (all 1.0 if UseEdgeLength=false).
#Returns NaN if V[1] is NaN, so that slot is left as it is.
if !isfinite(V[1])
    return NaN
end
total=0.0
totalL=0.0
for i=1:length(V)
    if isfinite(V[i])
        if !(isfinite(L[i]) && L[i]>0)
            error("Edge length $(L[i]) is not valid for an edge that has a value, do the Fe structs match this mesh?")
        end
        total+=V[i]*L[i]
        totalL+=L[i]
    end
end
return total/totalL
end

function EdgeLength(k,Idx,UseEdgeLength,FeP1P2S,FeP1P3S,FeP2P3S)
#Length of the free edge picked by GrabPointNew6: k=1 -> P2P3, k=2 -> P1P3, k=3 -> P1P2
if !UseEdgeLength
    return 1.0
end
if k==1
    return FeP2P3S.FeLe[Idx]
elseif k==2
    return FeP1P3S.FeLe[Idx]
elseif k==3
    return FeP1P2S.FeLe[Idx]
end
error("Empty row in InnerPoints")
end