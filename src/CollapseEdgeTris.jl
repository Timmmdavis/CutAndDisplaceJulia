function SortPointsOrderingAlongEdgeVector(FreePaPb,Pa,Pb,Pc)

#MdPnt to point Pb
FeM2Pb=Array{Float64,2}(undef, size(Pa));FeM2Pb=fill!(FeM2Pb, NaN)
Indx=findall(FreePaPb.FreeFlg)
#Check that Pb is always down the direction of FeEv - useful for cleaning tris later
FeM2Pb[Indx,:]=normr([Pb[Indx,1]-FreePaPb.FeMd[Indx,1] Pb[Indx,2]-FreePaPb.FeMd[Indx,2] Pb[Indx,3]-FreePaPb.FeMd[Indx,3]]);
for i=1:length(Indx)
    AllignFlag=dot(FreePaPb.FeEv[Indx[i],:],FeM2Pb[Indx[i],:]);
    if AllignFlag<0
        tmp=copy(Pb[Indx[i],:])
        Pb[Indx[i],:]=copy(Pa[Indx[i],:])
        Pa[Indx[i],:]=tmp
    end
end
return Pa,Pb,Pc
end


function CollapseEdgeTris(P1,P2,P3,MidPoint,FaceNormalVector)

#( Area,HalfPerimeter ) = CutAndDisplaceJulia.AreaOfTriangle3D( P1[:,1],P1[:,2],P1[:,3],P2[:,1],P2[:,2],P2[:,3],P3[:,1],P3[:,2],P3[:,3] );
#G=zeros(size(P1,1),1)
#G[114]=1; #172


###########################P1########################################
#RemoveAnyTrianglesWithMoreThan1EdgeFirst
#Number of edges:
(SortedTriangles,ConnectedEdge)=ConnectedConstraints(P1,P2,P3,MidPoint)
NoConnections=sum(SortedTriangles.!=0,dims=2)
Good=fill(true,size(NoConnections))
for i=1:length(NoConnections)
	if NoConnections[i]<2
		Good[i]=false
	end
end

P1new=zeros(sum(Good),3)
P2new=zeros(sum(Good),3)
P3new=zeros(sum(Good),3)
MidPointnew=zeros(sum(Good),3)
FaceNormalVectornew=zeros(sum(Good),3)
for i=1:sum(Good)
	if Good[i]==true
		P1new[i,:]=P1[i,:]
		P2new[i,:]=P2[i,:]
		P3new[i,:]=P3[i,:]
		MidPointnew[i,:]=MidPoint[i,:]
		FaceNormalVectornew[i,:]=FaceNormalVector[i,:]
	end
end
P1=copy(P1new)
P2=copy(P2new)
P3=copy(P3new)
MidPoint=copy(MidPointnew)
FaceNormalVector=copy(FaceNormalVectornew)


###########################P2########################################
#Sorting so the the last point described in FeP1P2S FeP1P3S FeP2P3S will always lie down the direction of FeEv
#(FeP1P2S,FeP1P3S,FeP2P3S)=GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector)
#(P1,P2,P3)=SortPointsOrderingAlongEdgeVector(FeP1P2S,P1,P2,P3);
#(P1,P3,P2)=SortPointsOrderingAlongEdgeVector(FeP1P3S,P1,P3,P2);
#(P2,P3,P1)=SortPointsOrderingAlongEdgeVector(FeP2P3S,P2,P3,P1);

(FeP1P2S,FeP1P3S,FeP2P3S)=GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector)
(SortedTriangles,ConnectedEdge)=ConnectedConstraints(P1,P2,P3,MidPoint)

###########################P3########################################
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


triindx=FirstTriIndx
OuterPointB=[NaN NaN NaN]
Collapse=false
removeIndx=NaN;
fillme=zeros(1,9)
newTris=fill(NaN,1,9)
EdgeTriInnerPoint=[0 0 0]
EdgeTriNonConnectedPoint=[0 0 0]
Next2EdgeTriNonConnectedPoint=[0 0 0]
Next2EdgeTriIndx=0;
Next2EdgeTriIndex=0;
RunOnceMore=false
for i=1:n_edges

	#@info triindx CurrentEdge 
	#@info P1[triindx,:] P2[triindx,:] P3[triindx,:]
	#@info CurrentPoint

	#Find Connected triangle in the correct direction
	for j=1:3

	    NeighbourIndex=SortedTriangles[triindx,j]
	    if NeighbourIndex==0
	        continue
	    end

	    #twoded
	    Next2EdgeTriIndx=NeighbourIndex;
	    
	    if EdgeTri[Next2EdgeTriIndx] #both edges and connected


	    	#Check that the next edge point is in the right direction
			(NewCurrentEdge)=LoopingRoundBoundaryKnownNextEdgeTri(triindx,CurrentEdge,CurrentPoint,TrailingPoint,Next2EdgeTriIndx,P1,P2,P3,FeP1P2S.FreeFlg,FeP1P3S.FreeFlg,FeP2P3S.FreeFlg)

			if NewCurrentEdge==0
				println("Not going the correct way")
				RunOnceMore=false
				continue #jump to next interation
			

			else #EdgePointsMatch

	    		(EdgeTriNonConnectedPoint,EdgeTriInnerPoint)=GetSpecficPointsCurrentTri(triindx,j,P1,P2,P3,FeP1P2S.FreeFlg,FeP1P3S.FreeFlg,FeP2P3S.FreeFlg)
	    		(Next2EdgeTriNonConnectedPoint)=GetSpecficPointsConnectedTri(triindx,j,Next2EdgeTriIndx,ConnectedEdge,P1,P2,P3,FeP1P2S.FreeFlg,FeP1P3S.FreeFlg,FeP2P3S.FreeFlg)

	    		removeIndx=[removeIndx triindx] #add to list of tris to remove

	    		OuterPointB=Next2EdgeTriNonConnectedPoint;

	    		println("Changing triindx to $Next2EdgeTriIndx")
	    		triindx=Next2EdgeTriIndx #update to the new tri

	    		CurrentEdge=NewCurrentEdge;

	    		#Flag to say we are collapsing triangles
	    		Collapse=true

	    		RunOnceMore=true
	    		break #drop out of loop

	    	end

	    else

	    	RunOnceMore=false

	    end

	end
	if RunOnceMore==true
		continue
	end

	if Collapse==true

		#If we get here we know that there is no edge connection to the current triangle
		#Adding to list of new triangles
		fillme[1:3].=EdgeTriInnerPoint;
		fillme[4:6].=Next2EdgeTriNonConnectedPoint;
		fillme[7:9].=EdgeTriNonConnectedPoint;
		#make sure point ordering matches that of the normal direction
		AvgVect=(FaceNormalVector[i,:].+FaceNormalVector[Next2EdgeTriIndx,:])./2;
		NewTriNormalVector = CreateTriangleNormal( fillme[1:3],fillme[4:6],fillme[7:9] );
		AllignFlag = dot(AvgVect,NewTriNormalVector);
		if AllignFlag<0
			#flip order
		    fillme[1:3].=EdgeTriNonConnectedPoint ;
		    fillme[7:9].=EdgeTriInnerPoint;
		end  

		newTris=[newTris; copy(fillme) ]
		Collapse=false

		#update the currentpoint
		TrailingPoint=CurrentPoint;
		CurrentPoint=Next2EdgeTriNonConnectedPoint;

	end

	#FindNextPointAlongFromThisTriangleEdgeAndSetNewTriIndex
	
	(triindx,CurrentPoint,TrailingPoint,CurrentEdge)=LoopingRoundBoundary(triindx,CurrentPoint,TrailingPoint,CurrentEdge,P1,P2,P3,FeP1P2S,FeP1P3S,FeP2P3S)
	Collapse=false

end

@info newTris
@info removeIndx

end


function GetSpecficPointsCurrentTri(triindx,j,P1,P2,P3,P1P2FreeFlg,P1P3FreeFlg,P2P3FreeFlg)
#Getting connected edges of each tri (P1 P2 or P3)
#j defines if its Col1=P1P2 Col2=P2P3 Col3=P1P3

if j==1
    #EdgeTriConnectedEdge=12
    EdgeTriNonConnectedPoint=P3[triindx,:]
    if P2P3FreeFlg[triindx]
        EdgeTriInnerPoint=P1[triindx,:]
    else
        EdgeTriInnerPoint=P2[triindx,:]
    end
elseif j==2
    #EdgeTriConnectedEdge=23
    EdgeTriNonConnectedPoint=P1[triindx,:]
    if P1P2FreeFlg[triindx]
        EdgeTriInnerPoint=P3[triindx,:]
    else
        EdgeTriInnerPoint=P2[triindx,:]
    end                    
elseif j==3
    #EdgeTriConnectedEdge=13
    EdgeTriNonConnectedPoint=P2[triindx,:]
    if P1P2FreeFlg[triindx]
        EdgeTriInnerPoint=P3[triindx,:]
    else
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
    #if P2P3FreeFlg[Next2EdgeTri]
    #    Next2EdgeTriInnerPoint=P1[Next2EdgeTri,:]
    #else
    #    Next2EdgeTriInnerPoint=P2[Next2EdgeTri,:]
    #end                    
elseif Next2EdgeTriEdge==23
    Next2EdgeTriNonConnectedPoint=P1[Next2EdgeTriIndx,:]
    #if P1P2FreeFlg[Next2EdgeTri]
    #    Next2EdgeTriInnerPoint=P3[Next2EdgeTri,:]
    #else
    #    Next2EdgeTriInnerPoint=P2[Next2EdgeTri,:]
    #end          
elseif Next2EdgeTriEdge==13
    Next2EdgeTriNonConnectedPoint=P2[Next2EdgeTriIndx,:]
    #if P1P2FreeFlg[Next2EdgeTri]
    #    Next2EdgeTriInnerPoint=P3[Next2EdgeTri,:]
    #else
    #    Next2EdgeTriInnerPoint=P1[Next2EdgeTri,:]
    #end   
end


return Next2EdgeTriNonConnectedPoint

end


function LoopingRoundBoundaryKnownNextEdgeTri(triindx,CurrentEdge,CurrentPoint,TrailingPoint,JoinedTriIndex,P1,P2,P3,P1P2FreeFlg,P1P3FreeFlg,P2P3FreeFlg)


	#reset
	NewCurrentEdge=0

	CurrentInP1=CurrentPoint==P1[JoinedTriIndex];
	if CurrentInP1==true & P1P2FreeFlg[JoinedTriIndex]==true
		#If the other point not the trailing point then go ahead
		if TrailingPoint!=P2[JoinedTriIndex]; 
			NewCurrentEdge=12
		end

	elseif CurrentInP1==true & P1P3FreeFlg[JoinedTriIndex]==true
		#If the other point not the trailing point then go ahead
		if TrailingPoint!=P3[JoinedTriIndex]; 
			NewCurrentEdge=13
		end		
	end

	InP2=CurrentPoint==P2[JoinedTriIndex];
	if InP2==true & P1P2FreeFlg[JoinedTriIndex]==true
		#If the other point not the trailing point then go ahead
		if TrailingPoint!=P1[JoinedTriIndex]; 
			NewCurrentEdge=12
		end

	elseif InP2==true & P2P3FreeFlg[JoinedTriIndex]==true
		#If the other point not the trailing point then go ahead
		if TrailingPoint!=P3[JoinedTriIndex]; 
			NewCurrentEdge=23
		end		

	end	

	InP3=CurrentPoint==P3[JoinedTriIndex];
	if InP3==true & P1P3FreeFlg[JoinedTriIndex]==true
		#If the other point not the trailing point then go ahead
		if TrailingPoint!=P1[JoinedTriIndex]; 
			NewCurrentEdge=13
		end			

	elseif InP3==true & P2P3FreeFlg[JoinedTriIndex]==true
		#If the other point not the trailing point then go ahead
		if TrailingPoint!=P2[JoinedTriIndex]; 
			NewCurrentEdge=23
		end					

	end	

return NewCurrentEdge

end

function LoopingRoundBoundary(triindx,CurrentPoint,TrailingPoint,CurrentEdge,P1,P2,P3,FeP1P2S,FeP1P3S,FeP2P3S)


	P1P2FreeFlg=FeP1P2S.FreeFlg;
	P1P3FreeFlg=FeP1P3S.FreeFlg;
	P2P3FreeFlg=FeP2P3S.FreeFlg;

	

	JoinedTriIndex=0
	InP1=ismember(P1,vec(CurrentPoint));
	if any(InP1) 
		Inside=findall(InP1)
		for i=1:length(Inside)


			if P1P2FreeFlg[Inside[i]]==true #Now check its a free edge
		
				(NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint)=
				GetNewEdgePoint(Inside[i],CurrentPoint,TrailingPoint,P1,P2,12)


			elseif P1P3FreeFlg[Inside[1]]==true

				(NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint)=
				GetNewEdgePoint(Inside[i],CurrentPoint,TrailingPoint,P1,P3,13)


			end

		end
	end


	InP2=ismember(P2,vec(CurrentPoint));
	if any(InP2) 
		Inside=findall(InP2)
		for i=1:length(Inside)

			if P1P2FreeFlg[Inside[i]]==true #Now check its a free edge

				(NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint)=
				GetNewEdgePoint(Inside[i],CurrentPoint,TrailingPoint,P1,P2,12)


			elseif P2P3FreeFlg[Inside[1]]==true

				(NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint)=
				GetNewEdgePoint(Inside[i],CurrentPoint,TrailingPoint,P2,P3,23)

			end

		end
	end
	

	InP3=ismember(P3,vec(CurrentPoint));
	if any(InP3) 
		Inside=findall(InP3)
		for i=1:length(Inside)



			if P1P3FreeFlg[Inside[i]]==true #Now check its a free edge


				(NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint)=
				GetNewEdgePoint(Inside[i],CurrentPoint,TrailingPoint,P1,P3,13)



			elseif P2P3FreeFlg[Inside[1]]==true

				(NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint)=
				GetNewEdgePoint(Inside[i],CurrentPoint,TrailingPoint,P2,P3,23)

			end

		end
	end	

	triindx=NewTriIndex;
	CurrentPoint=NewCurrentPoint;
	TrailingPoint=NewTrailingPoint;
	CurrentEdge=NewCurrentEdge;
	


	@info any(InP1) any(InP2) any(InP3) JoinedTriIndex CurrentEdge

	#V1=normr([P1[Inside[1],1]-CurrentPoint[1] P1[Inside[1],2]-CurrentPoint[2] P1[Inside[1],3]-CurrentPoint[3]])
				#AllignFlag=dot(V1,Ev)

return triindx,CurrentPoint,TrailingPoint,CurrentEdge

end


function GetNewEdgePoint(TestIndex,CurrentPoint,TrailingPoint,Pa,Pb,EdgeNo)
#EdgeNo=12;

println("Information ----------------->")
@info Pa[TestIndex,:]!=CurrentPoint Pb[TestIndex,:]!=CurrentPoint

if Pb[TestIndex,:]!=CurrentPoint && Pb[TestIndex,:]!=TrailingPoint

	NewTriIndex=TestIndex
	NewCurrentEdge=EdgeNo
	NewCurrentPoint=Pb[JoinedTriIndex,:]
	NewTrailingPoint=CurrentPoint

elseif Pa[TestIndex,:]!=CurrentPoint && Pa[TestIndex,:]!=TrailingPoint

	NewTriIndex=TestIndex
	NewCurrentEdge=EdgeNo
	NewCurrentPoint=Pa[JoinedTriIndex,:]
	NewTrailingPoint=CurrentPoint
end



return NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint

end