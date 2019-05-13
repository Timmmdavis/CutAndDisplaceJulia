function CollapseEdgeTris(P1,P2,P3,MidPoint,FaceNormalVector)

#=
#Remove any slither tris
( Area,HalfPerimeter ) = CutAndDisplaceJulia.AreaOfTriangle3D( P1[:,1],P1[:,2],P1[:,3],P2[:,1],P2[:,2],P2[:,3],P3[:,1],P3[:,2],P3[:,3] );
Good=vec(fill(true,length(Area)))
tol=mean(Area)/3
for i=1:length(Good)
	if Area[i]<tol
		PGood=false
	end
end
P1=copy(P1[Good,1:3])
P2=copy(P2[Good,1:3])
P3=copy(P3[Good,1:3])
(Points,Triangles)=CreateTrianglesPointsFromP1P2P3(P1,P2,P3)
(FaceNormalVector,MidPoint) = CreateFaceNormalAndMidPoint(Points,Triangles)
=#

## The aim is to collapse edge triangles that share an inner point
# We do this by looping round and checking if tris share a point inside the surface
#
#   T   C 	 	       T   C 		          T   ---> C
# ¯.¯\¯¯|¯¯/¯.¯  ¯.¯\¯¯|¯¯/¯.¯       ¯.¯\¯¯ ¯¯/¯.¯
# ....\ | /....  ....\ | /....  -- > ....\   /....
# .....\|/.....  .....\|/.....       .....\ /.....
# ......I......  ......I.......       .....I.......     
#	    
#  
# T=TrailingPoint
# C=CurrentPoint
# I=InnerPoint

(FeP1P2S,FeP1P3S,FeP2P3S)=GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector)
(SortedTriangles,ConnectedEdge)=ConnectedConstraints(P1,P2,P3,MidPoint)

#number of connected tris (Sorted tris rows not == to 0)
NoConnections=sum(SortedTriangles.!=0,dims=2)
EdgeTri=NoConnections.<3 #- edge tri

#Number of edges:
n_edges=sum([FeP1P2S.FreeFlg;FeP2P3S.FreeFlg;FeP1P3S.FreeFlg])

#Indx's of all the edges
FeP1P2Indxs=findall(FeP1P2S.FreeFlg);
FeP2P3Indxs=findall(FeP2P3S.FreeFlg);
FeP1P3Indxs=findall(FeP1P3S.FreeFlg);


if any(FeP1P2S.FreeFlg)
	#Start with First index in P1P2FreeFlg[:,1]
	CurrentEdge=12
	FirstTriIndx=FeP1P2Indxs[1];
	Direction=FeP1P2S.FeEv[FirstTriIndx] 
	EdgeMidPoint=FeP1P2S.FeMd[FirstTriIndx] 
	CurrentPoint=P2[FirstTriIndx,:]
	TrailingPoint=P1[FirstTriIndx,:]

elseif any(FeP2P3S.FreeFlg) #no edges in P1P2

	CurrentEdge=23
	FirstTriIndx=FeP2P3Indxs[1];
	Direction=FeP2P3S.FeEv[FirstTriIndx] 
	EdgeMidPoint=FeP2P3S.FeMd[FirstTriIndx]  
	CurrentPoint=P3[FirstTriIndx,:]
	TrailingPoint=P2[FirstTriIndx,:]

elseif any(FeP1P3S.FreeFlg)

	CurrentEdge=13
	FirstTriIndx=FeP1P3Indxs[1];
	Direction=FeP1P3S.FeEv[FirstTriIndx] 
	EdgeMidPoint=FeP1P3S.FeMd[FirstTriIndx]  
	CurrentPoint=P3[FirstTriIndx,:]
	TrailingPoint=P1[FirstTriIndx,:]

else
	error("It seems a little odd, your surface has no edges")
end

#The current index of the triangle on the outer edge
triindx=FirstTriIndx
#Flag to say we are inside part of loop that is collapsing triangles
Collapse=false
#Index's of triangles we will remove after loop
removeIndx=-1000; #To be removed
#Memoryless vector to fill with the list of new points
fillme=zeros(1,9)
#Appending to this vector everytime we add a new tri
newTris=fill(NaN,1,9)

EdgeTriInnerPoint=[0 0 0]
EdgeTriNonConnectedPoint=[0 0 0]
Next2EdgeTriNonConnectedPoint=[0 0 0]
Next2EdgeTriIndx=0;
Next2EdgeTriIndex=0;
RunOnceMore=false

for i=1:n_edges

	#Find Connected triangle in the correct direction
	for j=1:3

	    NeighbourIndex=SortedTriangles[triindx,j]
	    if NeighbourIndex==0
	        continue
	    end

	    #twoded
	    Next2EdgeTriIndx=NeighbourIndex;
	    
	    if EdgeTri[Next2EdgeTriIndx] #both edges and connected

	    	removeIndx=[removeIndx triindx] #add to list of tris to remove

	    	#Check that the next edge point is in the right direction
			(NewCurrentEdge)=
			LoopingRoundBoundaryKnownNextEdgeTri(triindx,CurrentEdge,CurrentPoint,TrailingPoint,Next2EdgeTriIndx,P1,P2,P3,
				FeP1P2S.FreeFlg,FeP1P3S.FreeFlg,FeP2P3S.FreeFlg)


			if NewCurrentEdge==0

				#Not going the correct way
				RunOnceMore=false

				continue #jump to next interation
			
			else #EdgePointsMatch


	    		(EdgeTriNonConnectedPointF,EdgeTriInnerPointF)=
	    		GetSpecficPointsCurrentTri(triindx,j,P1,P2,P3,FeP1P2S.FreeFlg,FeP1P3S.FreeFlg,FeP2P3S.FreeFlg)

	    		(Next2EdgeTriNonConnectedPoint)=
	    		GetSpecficPointsConnectedTri(triindx,j,Next2EdgeTriIndx,ConnectedEdge,P1,P2,P3,
	    			FeP1P2S.FreeFlg,FeP1P3S.FreeFlg,FeP2P3S.FreeFlg)

	    		#Update to the new tri we are working on
	    		triindx=Next2EdgeTriIndx 
	    		CurrentEdge=NewCurrentEdge;
	    		TrailingPoint=copy(CurrentPoint);
	    		CurrentPoint=copy(Next2EdgeTriNonConnectedPoint);

	    		if Collapse==false
	    			EdgeTriInnerPoint=EdgeTriInnerPointF
	    			EdgeTriNonConnectedPoint=EdgeTriNonConnectedPointF
	    		end
	    		
	    		Collapse=true 		#We are removing edges with shared inner Points
	    		RunOnceMore=true 	#We run the outer i loop once more before assuming the next tri is not a edge shared with this inner point
	    		break #drop out of inner j loop

	    	end

	    else

	    	RunOnceMore=false

	    end

	end

	#run outer i loop again
	if RunOnceMore==true
		continue
	end

	if Collapse==true

		#If we get here we know that there is no edge connection to the current triangle
		#Adding to list of new triangles
		fillme[1:3].=EdgeTriInnerPoint;
		fillme[4:6].=CurrentPoint;
		fillme[7:9].=EdgeTriNonConnectedPoint;

		
		#make sure point ordering matches that of the normal direction
		AvgVect=(FaceNormalVector[i,:].+FaceNormalVector[Next2EdgeTriIndx,:])./2;
		NewTriNormalVector = CreateTriangleNormal( fillme[1:3],fillme[4:6],fillme[7:9] );
		AllignFlag = dot(AvgVect,NewTriNormalVector);
		if AllignFlag<0
			#flip order
		    fillme[1:3].=EdgeTriNonConnectedPoint;
		    fillme[7:9].=EdgeTriInnerPoint;
		end  
		
		newTris=[newTris; copy(fillme) ]
		Collapse=false

	end

	#Find next outer triangle along from this triangle (they share CurrentPoint on the outer edge)
	(triindx,CurrentPoint,TrailingPoint,CurrentEdge)=
	LoopingRoundBoundary(triindx,CurrentPoint,TrailingPoint,CurrentEdge,P1,P2,P3,FeP1P2S,FeP1P3S,FeP2P3S)
	
end

return newTris,removeIndx

end


function GetSpecficPointsCurrentTri(triindx,j,P1,P2,P3,P1P2FreeFlg,P1P3FreeFlg,P2P3FreeFlg)
#Getting connected edges of each tri (P1 P2 or P3)
#j defines if its Col1=P1P2 Col2=P2P3 Col3=P1P3

if j==1
    #EdgeTriConnectedEdge=12
    EdgeTriNonConnectedPoint=P3[triindx,:]
    if P2P3FreeFlg[triindx]
        EdgeTriInnerPoint=P1[triindx,:]
    else #13free
        EdgeTriInnerPoint=P2[triindx,:]
    end
elseif j==2
    #EdgeTriConnectedEdge=23
    EdgeTriNonConnectedPoint=P1[triindx,:]
    if P1P2FreeFlg[triindx]
        EdgeTriInnerPoint=P3[triindx,:]
    else #13free
        EdgeTriInnerPoint=P2[triindx,:]
    end                    
elseif j==3
    #EdgeTriConnectedEdge=13
    EdgeTriNonConnectedPoint=P2[triindx,:]
    if P1P2FreeFlg[triindx]
        EdgeTriInnerPoint=P3[triindx,:]
    else #23free
        EdgeTriInnerPoint=P1[triindx,:]
    end                       
end

return EdgeTriNonConnectedPoint,EdgeTriInnerPoint

end


function GetSpecficPointsConnectedTri(triindx,j,Next2EdgeTriIndx,ConnectedEdge,P1,P2,P3,P1P2FreeFlg,P1P3FreeFlg,P2P3FreeFlg)
#Getting connected edges of each tri (P1 P2 or P3)
#j defines if its Col1=P1P2 Col2=P2P3 Col3=P1P3

Next2EdgeTriEdge=ConnectedEdge[triindx,j]
if Next2EdgeTriEdge==12
    Next2EdgeTriNonConnectedPoint=P3[Next2EdgeTriIndx,:]
elseif Next2EdgeTriEdge==23
    Next2EdgeTriNonConnectedPoint=P1[Next2EdgeTriIndx,:]
elseif Next2EdgeTriEdge==13
    Next2EdgeTriNonConnectedPoint=P2[Next2EdgeTriIndx,:]
end

return Next2EdgeTriNonConnectedPoint

end


function LoopingRoundBoundaryKnownNextEdgeTri(triindx,CurrentEdge,CurrentPoint,TrailingPoint,JoinedTriIndex,
	P1,P2,P3,P1P2FreeFlg,P1P3FreeFlg,P2P3FreeFlg)


	#reset
	NewCurrentEdge=0

	CurrentInP1=CurrentPoint==P1[JoinedTriIndex,:];
	if CurrentInP1==true & P1P2FreeFlg[JoinedTriIndex]==true
		#If the other point not the trailing point then go ahead
		if TrailingPoint!=P2[JoinedTriIndex,:]; 
			NewCurrentEdge=12
		end

	elseif CurrentInP1==true & P1P3FreeFlg[JoinedTriIndex]==true
		#If the other point not the trailing point then go ahead
		if TrailingPoint!=P3[JoinedTriIndex,:]; 
			NewCurrentEdge=13
		end		
	end

	CurrentInP2=CurrentPoint==P2[JoinedTriIndex,:];
	if CurrentInP2==true & P1P2FreeFlg[JoinedTriIndex]==true
		#If the other point not the trailing point then go ahead
		if TrailingPoint!=P1[JoinedTriIndex,:]; 
			NewCurrentEdge=12
		end

	elseif CurrentInP2==true & P2P3FreeFlg[JoinedTriIndex]==true
		#If the other point not the trailing point then go ahead
		if TrailingPoint!=P3[JoinedTriIndex,:]; 
			NewCurrentEdge=23
		end		

	end	

	CurrentInP3=CurrentPoint==P3[JoinedTriIndex,:];
	if CurrentInP3==true & P1P3FreeFlg[JoinedTriIndex]==true
		#If the other point not the trailing point then go ahead
		if TrailingPoint!=P1[JoinedTriIndex,:]; 
			NewCurrentEdge=13
		end			

	elseif CurrentInP3==true & P2P3FreeFlg[JoinedTriIndex]==true
		#If the other point not the trailing point then go ahead
		if TrailingPoint!=P2[JoinedTriIndex,:]; 
			NewCurrentEdge=23
		end					

	end	

return NewCurrentEdge

end

function LoopingRoundBoundary(triindx,CurrentPoint,TrailingPoint,CurrentEdge,
								P1,P2,P3,FeP1P2S,FeP1P3S,FeP2P3S)


	P1P2FreeFlg=FeP1P2S.FreeFlg;
	P1P3FreeFlg=FeP1P3S.FreeFlg;
	P2P3FreeFlg=FeP2P3S.FreeFlg;


	
	NewTriIndex=0 #NewTriIndex
	NewCurrentEdge=0 #NewCurrentEdge
	NewCurrentPoint=[0 0 0] #NewCurrentPoint
	NewTrailingPoint=[0 0 0] #NewTrailingPoint


	InP1=ismember(P1,vec(CurrentPoint));
	if any(InP1) 
		Inside=findall(InP1)
		for i=1:length(Inside)


			if P1P2FreeFlg[Inside[i]]==true #Now check its a free edge
		
				(NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint)=
				GetNewEdgePoint(Inside[i],CurrentPoint,TrailingPoint,P1,P2,12,
					NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint)

			end
			if P1P3FreeFlg[Inside[i]]==true

				(NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint)=
				GetNewEdgePoint(Inside[i],CurrentPoint,TrailingPoint,P1,P3,13,
					NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint)

			end

		end
	end


	InP2=ismember(P2,vec(CurrentPoint));
	if any(InP2) 
		Inside=findall(InP2)
		println(Inside)

		for i=1:length(Inside)


			if P1P2FreeFlg[Inside[i]]==true #Now check its a free edge

				(NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint)=
				GetNewEdgePoint(Inside[i],CurrentPoint,TrailingPoint,P1,P2,12,
					NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint)

			end
			if P2P3FreeFlg[Inside[i]]==true

				(NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint)=
				GetNewEdgePoint(Inside[i],CurrentPoint,TrailingPoint,P2,P3,23,
					NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint)

			end

		end
	end


	InP3=ismember(P3,vec(CurrentPoint));
	if any(InP3) 
		Inside=findall(InP3)
		for i=1:length(Inside)

			if P1P3FreeFlg[Inside[i]]==true #Now check its a free edge

				(NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint)=
				GetNewEdgePoint(Inside[i],CurrentPoint,TrailingPoint,P1,P3,13,
					NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint)

			end
			if P2P3FreeFlg[Inside[i]]==true

				(NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint)=
				GetNewEdgePoint(Inside[i],CurrentPoint,TrailingPoint,P2,P3,23,
					NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint)

			end

		end
	end	

	triindx=NewTriIndex;
	CurrentPoint=NewCurrentPoint;
	TrailingPoint=NewTrailingPoint;
	CurrentEdge=NewCurrentEdge;
	

return triindx,CurrentPoint,TrailingPoint,CurrentEdge

end


function GetNewEdgePoint(TestIndex,CurrentPoint,TrailingPoint,Pa,Pb,EdgeNo,NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint)

if Pb[TestIndex,:]!=CurrentPoint && Pb[TestIndex,:]!=TrailingPoint

	NewTriIndex=TestIndex
	NewCurrentEdge=EdgeNo
	NewCurrentPoint=Pb[TestIndex,:]
	NewTrailingPoint=CurrentPoint

elseif Pa[TestIndex,:]!=CurrentPoint && Pa[TestIndex,:]!=TrailingPoint

	NewTriIndex=TestIndex
	NewCurrentEdge=EdgeNo
	NewCurrentPoint=Pa[TestIndex,:]
	NewTrailingPoint=CurrentPoint
end

return NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint

end