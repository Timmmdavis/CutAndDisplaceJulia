function MovingAverageOfStressIntensity(avgeverynth,P1,P2,P3,FaceNormalVector,MidPoint,FeP1P2S,FeP1P3S,FeP2P3S)

#avgeverynth - if 3 averages from the surrounding triangles
#		 if 5 averages from two tris either side


for z=1:length(Variables)
	Variable=Variables[z];
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
		    	Vtrailing=FeP2P3S.StrainEnergy[Idx_trailing]
		    elseif k==2
		    	Vtrailing=FeP1P3S.StrainEnergy[Idx_trailing]
		    elseif k==3
		    	Vtrailing=FeP1P2S.StrainEnergy[Idx_trailing]
		    end    	
		    (k,Idx_future) =GrabPointNew6(InnerPoints,P1,P2,P3,j_future)
			if k==1
		    	Vfuture=FeP2P3S.StrainEnergy[Idx_future]
		    elseif k==2
		    	Vfuture=FeP1P3S.StrainEnergy[Idx_future]
		    elseif k==3
		    	Vfuture=FeP1P2S.StrainEnergy[Idx_future]
		    end   
		    (k,Idx_current) =GrabPointNew6(InnerPoints,P1,P2,P3,j_current)
			if k==1
				avg=(Vtrailing+FeP2P3S.StrainEnergy[Idx_current]+Vfuture)/3;
				if isnan(avg)
					continue
				end
				FeP2P3S.StrainEnergy[Idx_current]=avg
		    elseif k==2
		    	avg=(Vtrailing+FeP1P3S.StrainEnergy[Idx_current]+Vfuture)/3
				if isnan(avg)
					continue
				end
		    	FeP1P3S.StrainEnergy[Idx_current]=(Vtrailing+FeP1P3S.StrainEnergy[Idx_current]+Vfuture)/3
		    elseif k==3
		    	avg=(Vtrailing+FeP1P2S.StrainEnergy[Idx_current]+Vfuture)/3    	
		    	if isnan(avg)
					continue
				end
		    	FeP1P2S.StrainEnergy[Idx_current]=(Vtrailing+FeP1P2S.StrainEnergy[Idx_current]+Vfuture)/3    	
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
		    	Vtrailing=FeP2P3S.StrainEnergy[Idx_trailing]
		    elseif k==2
		    	Vtrailing=FeP1P3S.StrainEnergy[Idx_trailing]
		    elseif k==3
		    	Vtrailing=FeP1P2S.StrainEnergy[Idx_trailing]
		    end   
			(k,Idx_trailing) =GrabPointNew6(InnerPoints,P1,P2,P3,j_trailing2)
		    if k==1
		    	Vtrailing2=FeP2P3S.StrainEnergy[Idx_trailing]
		    elseif k==2
		    	Vtrailing2=FeP1P3S.StrainEnergy[Idx_trailing]
		    elseif k==3
		    	Vtrailing2=FeP1P2S.StrainEnergy[Idx_trailing]
		    end   
		    (k,Idx_future) =GrabPointNew6(InnerPoints,P1,P2,P3,j_future)
			if k==1
		    	Vfuture=FeP2P3S.StrainEnergy[Idx_future]
		    elseif k==2
		    	Vfuture=FeP1P3S.StrainEnergy[Idx_future]
		    elseif k==3
		    	Vfuture=FeP1P2S.StrainEnergy[Idx_future]
		    end   
		    (k,Idx_future) =GrabPointNew6(InnerPoints,P1,P2,P3,j_future2)
			if k==1
		    	Vfuture2=FeP2P3S.StrainEnergy[Idx_future]
		    elseif k==2
		    	Vfuture2=FeP1P3S.StrainEnergy[Idx_future]
		    elseif k==3
		    	Vfuture2=FeP1P2S.StrainEnergy[Idx_future]
		    end   	    
		    (k,Idx_current) =GrabPointNew6(InnerPoints,P1,P2,P3,j_current)
			if k==1
				avg=(Vtrailing2+Vtrailing+FeP2P3S.StrainEnergy[Idx_current]+Vfuture+Vfuture2)/5;
				if isnan(avg)
					continue
				end
				FeP2P3S.StrainEnergy[Idx_current]=avg
		    elseif k==2
		    	avg=(Vtrailing2+Vtrailing+FeP1P3S.StrainEnergy[Idx_current]+Vfuture+Vfuture2)/5
				if isnan(avg)
					continue
				end
		    	FeP1P3S.StrainEnergy[Idx_current]=avg
		    elseif k==3
		    	avg=(Vtrailing2+Vtrailing+FeP1P2S.StrainEnergy[Idx_current]+Vfuture+Vfuture2)/5   	
		    	if isnan(avg)
					continue
				end
		    	FeP1P2S.StrainEnergy[Idx_current]=avg
		    end       
		    

		    #@info Idx_trailing Idx_current Idx_future
		    end
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