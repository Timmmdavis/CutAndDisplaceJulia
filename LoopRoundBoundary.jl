function LoopRoundBoundary(MidPoint,P1,P2,P3,Triangles,FaceNormalVector)

(P1P2FreeFlg,P2P3FreeFlg,P1P3FreeFlg)=EdgeConstraints(P1,P2,P3,MidPoint);
(SortedTriangles,ConnectedEdge)=ConnectedConstraints(P1,P2,P3,MidPoint)
(FeP1P2S,FeP1P3S,FeP2P3S)=GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector)

#Number of edges:
n_edges=sum([FeP1P2S.FreeFlag;FeP2P3.FreeFlag;FeP1P3.FreeFlag])

#Indx's of all the edges
FeP1P2Indxs=findall(FeP1P2S.FreeFlag);
FeP2P3Indxs=findall(FeP2P3S.FreeFlag);
FeP1P3Indxs=findall(FeP1P3S.FreeFlag);

if any(FeP1P2S.FreeFlag)
	#Start with First index in P1P2FreeFlg[:,1]
	
	FirstTriIndx=FeP1P2Indxs[1];
	Direction=FeP1P2S.FeEvp[FirstTriIndx] 
	EdgeMidPoint=FeP1P2S.FeMd[FirstTriIndx] 

elseif any(FeP2P3.FreeFlag) #no edges in P1P2

	Indx=findall(FeP2P3.FreeFlag);
	FirstTriIndx=Indx[1];
	Direction=FeP2P3.FeEv[FirstTriIndx] 
	EdgeMidPoint=FeP2P3.FeMd[FirstTriIndx]  

elseif any(FeP1P3.FreeFlag)

	Indx=findall(FeP1P3.FreeFlag);
	FirstTriIndx=Indx[1];
	Direction=FeP1P3.FeEv[FirstTriIndx] 
	EdgeMidPoint=FeP1P3.FeMd[FirstTriIndx]  

else
	error("It seems a little odd, your surface has no edges")
end

NoConnections=sum(SortedTriangles.!=0,dims=2)
EdgeTri=NoConnections.<3 #- edge tri




triindx=FirstTriIndx
for i=1:n_edges


	(EdgeTriInnerPoint,
		EdgeTriNonConnectedPoint,
		Next2EdgeTriNonConnectedPoint,
		MultiShareInnerPoint,
		triindx,
		mergeme,
		removeindx)=
	TestEdgeConnectionsOfCurrentTri(
		ConnectedTris,
		triindx,
		EdgeTri,
		ConnectedEdge,
		EdgeMidPoint,
		Direction,
		P1,
		P2,
		P3,
		P1P2FreeFlg,
		P2P3FreeFlg,
		P1P3FreeFlg,
		MultiShareInnerPoint,
		FeP1P2Indxs,
		FeP2P3Indxs
		FeP1P3Indxs,
		EdgeTriInnerPoint,
		EdgeTriNonConnectedPoint,
		Next2EdgeTriNonConnectedPoint,
		removeindx)

	#=
	if mergeme==true

		#Adding to list of new triangles
		fillme[1:3].=EdgeTriInnerPoint;
		fillme[4:6].=Next2EdgeTriNonConnectedPoint;
		fillme[7:9].=EdgeTriNonConnectedPoint;
		#make sure point ordering matches that of the normal direction
		AvgVect=(FaceNormalVector[i,:].+FaceNormalVector[Next2EdgeTri,:])./2;
		NewTriNormalVector = CreateTriangleNormal( fillme[1:3],fillme[4:6],fillme[7:9] );
		C = dot(AvgVect,NewTriNormalVector);
		if C<0
		    fillme[1:3].=EdgeTriNonConnectedPoint ;
		    fillme[7:9].=EdgeTriInnerPoint;
		end  
		if append==0
		    newTris=copy(fillme)
		else
		    newTris=[newTris; copy(fillme) ]
		end

		#Adding to list of triangles to remove
		fillme2[1]=i;
		fillme2[2]=Next2EdgeTri;
		if append==0
		    removeIndx=fillme2
		else
		    removeIndx=[removeIndx; fillme2 ]
		end    

		#once we hit this the first time we want to add to the vectors
		append=1; 
	end
	=#

end
return removeindx
end


function TestEdgeConnectionsOfCurrentTri(
	ConnectedTris,
	triindx,
	EdgeTri,
	ConnectedEdge,
	EdgeMidPoint,
	Direction,
	P1,
	P2,
	P3,
	P1P2FreeFlg,
	P2P3FreeFlg,
	P1P3FreeFlg,
	MultiShareInnerPoint,
	FeP1P2Indxs,
	FeP2P3Indxs
	FeP1P3Indxs,
	EdgeTriInnerPoint,
	EdgeTriNonConnectedPoint,
	Next2EdgeTriNonConnectedPoint,
	removeindx)

#Find Connected triangle in the correct direction
for j=1:3

    Next2EdgeTri=ConnectedTris[triindx,j]
    if Next2EdgeTri==0
        continue
    end

    #add to the old remove indx before going forward
    removeindx=[removeindx,triindx];
    
    if EdgeTri[Next2EdgeTri] #both edges and connected

    	
    	

	    #Now find it its in the right direction

	    #The point numbering of the connected tris edge
	    ConnectedTriEdgeNumbering=ConnectedEdge[triindx,j]
	    if ConnectedTriEdgeNumbering==12

	    	NewEdgeMidPoint=FeP1P2S.FeMd[Next2EdgeTri]
	    	Next2EdgeTriNonConnectedPoint=P3[Next2EdgeTri,:]  

	    	if P2P3FreeFlg[Next2EdgeTri]
	    	    Next2EdgeTriInnerPoint=P1[Next2EdgeTri,:]
	    	    Next2EdgeTriConnectedPoint=P2[Next2EdgeTri,:]
	    	else
	    	    Next2EdgeTriInnerPoint=P2[Next2EdgeTri,:]
	    	    Next2EdgeTriConnectedPoint=P1[Next2EdgeTri,:]
	    	end     

	    elseif ConnectedTriEdgeNumbering==23

	    	NewEdgeMidPoint=FeP2P3.FeMd[Next2EdgeTri]
	    	Next2EdgeTriNonConnectedPoint=P1[Next2EdgeTri,:]  

	    	if P1P2FreeFlg[Next2EdgeTri]
	    	    Next2EdgeTriInnerPoint=P3[Next2EdgeTri,:]
	    	    Next2EdgeTriConnectedPoint=P2[Next2EdgeTri,:]
	    	else
	    	    Next2EdgeTriInnerPoint=P2[Next2EdgeTri,:]
	    	    Next2EdgeTriConnectedPoint=P3[Next2EdgeTri,:]
	    	end   

	    elseif ConnectedTriEdgeNumbering==13

	    	NewEdgeMidPoint=FeP1P3.FeMd[Next2EdgeTri]
	    	Next2EdgeTriNonConnectedPoint=P2[Next2EdgeTri,:]

	    	if P1P2FreeFlg[Next2EdgeTri]
	    	    Next2EdgeTriInnerPoint=P3[Next2EdgeTri,:]
	    	    Next2EdgeTriConnectedPoint=P1[Next2EdgeTri,:]
	    	else
	    	    Next2EdgeTriInnerPoint=P1[Next2EdgeTri,:]
	    	    Next2EdgeTriConnectedPoint=P3[Next2EdgeTri,:]
	    	end   

	    end

	    #Compute the distance between the two edge midpoints
	    DistanceBetweenMidPoints=sqrt((EdgeMidPoint[1]+NewEdgeMidPoint[1])^2+
	    							  (EdgeMidPoint[2]+NewEdgeMidPoint[2])^2+
	    							  (EdgeMidPoint[3]+NewEdgeMidPoint[3])^2)

	    #Move the first edge along our loop direction
	    EdgeMidPointDir=dot(EdgeMidPoint,Direction.*1e-10);

	    #check the distance again
	    Distance2=sqrt((EdgeMidPointDir[1]+NewEdgeMidPoint[1])^2+
	    			   (EdgeMidPointDir[2]+NewEdgeMidPoint[2])^2+
	    			   (EdgeMidPointDir[3]+NewEdgeMidPoint[3])^2)	   

		if Distance2>DistanceBetweenMidPoints						   
			println("Not the right edge triangle, wrong direction")
			continue #skip to next one
		end

		#Only do if multiple tris share inner point is off (if on these we want to retain
		#the first tri in this loop direction that had this inner point)
		if MultiShareInnerPoint=false 

			#Getting connected edges of each tri (P1 P2 or P3)
			#j defines if its Col1=P1P2 Col2=P2P3 Col3=P1P3
			if j==1
			    #EdgeTriEdge=12
			    EdgeTriNonConnectedPoint=P3[triindx,:]
			    if P2P3FreeFlg[triindx]
			        EdgeTriInnerPoint=P1[triindx,:]
			    else
			        EdgeTriInnerPoint=P2[triindx,:]
			    end
			elseif j==2
			    #EdgeTriEdge=23
			    EdgeTriNonConnectedPoint=P1[triindx,:]
			    if P1P2FreeFlg[triindx]
			        EdgeTriInnerPoint=P3[triindx,:]
			    else
			        EdgeTriInnerPoint=P2[triindx,:]
			    end                    
			elseif j==3
			    #EdgeTriEdge=13
			    EdgeTriNonConnectedPoint=P2[triindx,:]
			    if P1P2FreeFlg[triindx]
			        EdgeTriInnerPoint=P3[triindx,:]
			    else
			        EdgeTriInnerPoint=P1[triindx,:]
			    end                       
			end

			if all(EdgeTriInnerPoint.==Next2EdgeTriInnerPoint)
				MultiShareInnerPoint=true
			end

		end

		#Set index to the next triangle along
		CurrentTri=Next2EdgeTri;


	else

		#Its not also an edge tri!!!!
		MultiShareInnerPoint=false
		
		#Find the next index and return this!!!
		Next2EdgeTriConnectedPoint
		ConnectedIndx=ismember(P1[FeP1P2Indxs,:],Next2EdgeTriConnectedPoint)
		if any(ConnectedIndx)
			Loc=findall(ConnectedIndx)
			CurrentTri=FeP1P2Indxs[Loc]
			continue
		end
		ConnectedIndx=ismember(P2[FeP1P2Indxs,:],Next2EdgeTriConnectedPoint)
		if any(ConnectedIndx)
			Loc=findall(ConnectedIndx)
			CurrentTri=FeP1P2Indxs[Loc]
			continue
		end
		ConnectedIndx=ismember(P2[FeP2P3Indxs,:],Next2EdgeTriConnectedPoint)
		if any(ConnectedIndx)
			Loc=findall(ConnectedIndx)
			CurrentTri=FeP2P3Indxs[Loc]
			continue
		end		
		ConnectedIndx=ismember(P3[FeP2P3Indxs,:],Next2EdgeTriConnectedPoint)
		if any(ConnectedIndx)
			Loc=findall(ConnectedIndx)
			CurrentTri=FeP2P3Indxs[Loc]
			continue
		end				
		ConnectedIndx=ismember(P1[FeP1P3Indxs,:],Next2EdgeTriConnectedPoint)
		if any(ConnectedIndx)
			Loc=findall(ConnectedIndx)
			CurrentTri=FeP1P3Indxs[Loc]
			continue
		end					
		ConnectedIndx=ismember(P3[FeP1P3Indxs,:],Next2EdgeTriConnectedPoint)
		if any(ConnectedIndx)
			Loc=findall(ConnectedIndx)
			CurrentTri=FeP1P3Indxs[Loc]
			continue
		end			

	end

end

return EdgeTriInnerPoint,EdgeTriNonConnectedPoint,Next2EdgeTriNonConnectedPoint,MultiShareInnerPoint,CurrentTri,removeindx


end