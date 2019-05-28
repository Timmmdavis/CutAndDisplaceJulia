function GetSortedEdgesOfMeshList(P1,P2,P3,FaceNormalVector,MidPoint)
#Returns some vectors:

#UniqueEdges
#Step ranges of the edge in the lists below

#LeadingPoints,TrailingPoints,InnerPoints
#n*3 arrays where col 1 is P1 col2 P2 etc. Rows not equal to 0 are the index's

( Area,HalfPerimeter ) = CutAndDisplaceJulia.AreaOfTriangle3D( P1[:,1],P1[:,2],P1[:,3],P2[:,1],P2[:,2],P2[:,3],P3[:,1],P3[:,2],P3[:,3] );
(FeP1P2S,FeP1P3S,FeP2P3S)=GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector)
(SortedTriangles,ConnectedEdge)=ConnectedConstraints(P1,P2,P3,MidPoint)
#number of connected tris (Sorted tris rows not == to 0)
NoConnections=sum(SortedTriangles.!=0,dims=2)
EdgeTri=NoConnections.<3 #- edge tri
TotalNoFreeEdges=sum(SortedTriangles.==0)

#Number of edges:
n_edges=sum([FeP1P2S.FreeFlg;FeP2P3S.FreeFlg;FeP1P3S.FreeFlg])

#Vector where going down through the rows we loop through the boundary with the index's of P1, P2 and P3
#[P1  P2 P3  ]
#[NaN 678 NaN]
#[NaN NaN 50 ]
#[NaN 5   NaN]
LeadingPoints=	fill(0,TotalNoFreeEdges,3)
TrailingPoints=	fill(0,TotalNoFreeEdges,3)
InnerPoints=	fill(0,TotalNoFreeEdges,3)


Counter=0
UniqueEdges=[]
#Pre alloc
CurrentPoint=[NaN NaN NaN]# reset
FirstPnt=[NaN NaN NaN]# reset
TrailingPoint=[NaN NaN NaN]# reset
InnerPoint=[NaN NaN NaN]# reset
triindx=NaN
CurrentEdge=NaN

#While sorted triangles still has edges we have not passed
while sum(SortedTriangles.!=0)!=length(SortedTriangles) 

	#Getting the trailing and leading pnt - finding the first triangle that is an edge in the list
	#We change edges we have passed in the list to -1 leaving only edges we have not seen in the list as '0'
	breaki=false
	for i=1:size(SortedTriangles,1) #for each tri
		for j=1:size(SortedTriangles,2) #for each potential connection

			if SortedTriangles[i,j]==0 #Its a tri with a missing connection (free edge)

				CurrentPoint=[NaN NaN NaN]# reset
				triindx=i; #Current tri with free edge 
				if j==1

					CurrentEdge=12
					Direction=FeP1P2S.FeEv[triindx] 
					EdgeMidPoint=FeP1P2S.FeMd[triindx] 
					FirstPnt=P2[triindx,:]
					TrailingPoint=P1[triindx,:]

				elseif j==2

					CurrentEdge=23
					Direction=FeP2P3S.FeEv[triindx] 
					EdgeMidPoint=FeP2P3S.FeMd[triindx]  
					FirstPnt=P3[triindx,:]
					TrailingPoint=P2[triindx,:]

				elseif j==3

					CurrentEdge=13
					Direction=FeP1P3S.FeEv[triindx] 
					EdgeMidPoint=FeP1P3S.FeMd[triindx]  
					FirstPnt=P3[triindx,:]
					TrailingPoint=P1[triindx,:]

				else
					error("It seems a little odd, your surface has no edges")
				end

				breaki=true
				break
			end
		end
		#Leave i loop too
		if breaki==true
			break
		end
	end

	#Loop round this boundary until we rehit the first point 
	while CurrentPoint!=FirstPnt

		if isnan(CurrentPoint[1])
			CurrentPoint=FirstPnt
		end

		#Find next outer triangle along from this triangle (they share CurrentPoint on the outer edge)
		(triindx,CurrentPoint,TrailingPoint,CurrentEdge,InnerPoint,FrontPointLoc,BackPointLoc,InnerPointLoc)=
		LoopingRoundBoundaries(triindx,CurrentPoint,TrailingPoint,CurrentEdge,P1,P2,P3,FeP1P2S,FeP1P3S,FeP2P3S)

		#FrontPointLoc BackPointLoc InnerPointLoc - In the current triindx which
		#points (from P1P2P3) represent the front, back and inner point
		Counter+=1; #Indx inside our lists : TrailingPoints LeadingPoints etc

		if Counter>TotalNoFreeEdges
			error("Irk")
		end

		if FrontPointLoc==BackPointLoc || FrontPointLoc==InnerPointLoc || InnerPointLoc==BackPointLoc
			error("Whats happening here")
		end

		if FrontPointLoc==1
			LeadingPoints[Counter,1]=triindx
		elseif FrontPointLoc==2
			LeadingPoints[Counter,2]=triindx
		elseif FrontPointLoc==3
			LeadingPoints[Counter,3]=triindx
		end
		if BackPointLoc==1
			TrailingPoints[Counter,1]=triindx
		elseif BackPointLoc==2
			TrailingPoints[Counter,2]=triindx
		elseif BackPointLoc==3
			TrailingPoints[Counter,3]=triindx
		end
		if InnerPointLoc==1
			InnerPoints[Counter,1]=triindx
		elseif InnerPointLoc==2
			InnerPoints[Counter,2]=triindx
		elseif InnerPointLoc==3
			InnerPoints[Counter,3]=triindx
		end

		if CurrentEdge==12
			SortedTriangles[triindx,1]=-1
		elseif CurrentEdge==23
			SortedTriangles[triindx,2]=-1
		elseif CurrentEdge==13
			SortedTriangles[triindx,3]=-1
		end



	end

	if isempty(UniqueEdges)
		UniqueEdges=[1:Counter]
	else
		UniqueEdges=[UniqueEdges;[maximum(UniqueEdges[end])+1:Counter]]
	end
	CurrentPoint=[NaN NaN NaN] #reset


end

return UniqueEdges,LeadingPoints,TrailingPoints,InnerPoints
end



function LoopingRoundBoundaries(triindx,CurrentPoint,TrailingPoint,CurrentEdge,
								P1,P2,P3,FeP1P2S,FeP1P3S,FeP2P3S)


	P1P2FreeFlg=FeP1P2S.FreeFlg;
	P1P3FreeFlg=FeP1P3S.FreeFlg;
	P2P3FreeFlg=FeP2P3S.FreeFlg;


	
	NewTriIndex=0 #NewTriIndex
	NewCurrentEdge=0 #NewCurrentEdge
	NewCurrentPoint=[0. 0. 0.] #NewCurrentPoint
	NewTrailingPoint=[0. 0. 0.] #NewTrailingPoint
	NewInnerPoint=[0. 0. 0.]
	PaOrPb=0
	FrontPointLoc=0
	BackPointLoc=0
	InnerPointLoc=0

	InP1=ismember(P1,vec(CurrentPoint));
	if any(InP1) 
		Inside=findall(InP1)
		for i=1:length(Inside)


			if P1P2FreeFlg[Inside[i]]==true #Now check its a free edge
		
				(NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint,PaOrPb)=
				GetNewEdgePoint(Inside[i],CurrentPoint,TrailingPoint,P1,P2,12,
					NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint,PaOrPb)
				NewInnerPoint=P3[Inside[i],:]
				if PaOrPb==1
					FrontPointLoc=1;
					BackPointLoc=2;
					InnerPointLoc=3;
				elseif PaOrPb==2
					FrontPointLoc=2;
					BackPointLoc=1;
					InnerPointLoc=3;
				end

			end
			if P1P3FreeFlg[Inside[i]]==true

				(NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint,PaOrPb)=
				GetNewEdgePoint(Inside[i],CurrentPoint,TrailingPoint,P1,P3,13,
					NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint,PaOrPb)
				NewInnerPoint=P2[Inside[i],:]
				if PaOrPb==1
					FrontPointLoc=1;
					BackPointLoc=3;
					InnerPointLoc=2;
				elseif PaOrPb==2
					FrontPointLoc=3;
					BackPointLoc=1;
					InnerPointLoc=2;
				end

			end

		end
	end


	InP2=ismember(P2,vec(CurrentPoint));
	if any(InP2) 
		Inside=findall(InP2)

		for i=1:length(Inside)


			if P1P2FreeFlg[Inside[i]]==true #Now check its a free edge

				(NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint,PaOrPb)=
				GetNewEdgePoint(Inside[i],CurrentPoint,TrailingPoint,P1,P2,12,
					NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint,PaOrPb)
				NewInnerPoint=P3[Inside[i],:]
				if PaOrPb==1
					FrontPointLoc=1;
					BackPointLoc=2;
					InnerPointLoc=3;
				elseif PaOrPb==2
					FrontPointLoc=2;
					BackPointLoc=1;
					InnerPointLoc=3;
				end

			end
			if P2P3FreeFlg[Inside[i]]==true

				(NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint,PaOrPb)=
				GetNewEdgePoint(Inside[i],CurrentPoint,TrailingPoint,P2,P3,23,
					NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint,PaOrPb)
				NewInnerPoint=P1[Inside[i],:]
				if PaOrPb==1
					FrontPointLoc=2;
					BackPointLoc=3;
					InnerPointLoc=1;
				elseif PaOrPb==2
					FrontPointLoc=3;
					BackPointLoc=2;
					InnerPointLoc=1;
				end

			end

		end
	end


	InP3=ismember(P3,vec(CurrentPoint));
	if any(InP3) 
		Inside=findall(InP3)
		for i=1:length(Inside)

			if P1P3FreeFlg[Inside[i]]==true #Now check its a free edge

				(NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint,PaOrPb)=
				GetNewEdgePoint(Inside[i],CurrentPoint,TrailingPoint,P1,P3,13,
					NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint,PaOrPb)
				NewInnerPoint=P2[Inside[i],:]
				if PaOrPb==1
					FrontPointLoc=1;
					BackPointLoc=3;
					InnerPointLoc=2;
				elseif PaOrPb==2
					FrontPointLoc=3;
					BackPointLoc=1;
					InnerPointLoc=2;
				end

			end
			if P2P3FreeFlg[Inside[i]]==true

				(NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint,PaOrPb)=
				GetNewEdgePoint(Inside[i],CurrentPoint,TrailingPoint,P2,P3,23,
					NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint,PaOrPb)
				NewInnerPoint=P1[Inside[i],:]
				if PaOrPb==1
					FrontPointLoc=2;
					BackPointLoc=3;
					InnerPointLoc=1;
				elseif PaOrPb==2
					FrontPointLoc=3;
					BackPointLoc=2;
					InnerPointLoc=1;
				end

			end

		end
	end	

	triindx=NewTriIndex;
	CurrentPoint=NewCurrentPoint;
	TrailingPoint=NewTrailingPoint;
	CurrentEdge=NewCurrentEdge;
	InnerPoint=NewInnerPoint

return triindx,CurrentPoint,TrailingPoint,CurrentEdge,InnerPoint,FrontPointLoc,BackPointLoc,InnerPointLoc

end


function GetNewEdgePoint(TestIndex,CurrentPoint,TrailingPoint,Pa,Pb,EdgeNo,NewTriIndex,
	NewCurrentEdge,NewCurrentPoint,NewTrailingPoint,PaOrPb)

if Pb[TestIndex,:]!=CurrentPoint && Pb[TestIndex,:]!=TrailingPoint

	NewTriIndex=TestIndex
	NewCurrentEdge=EdgeNo
	NewCurrentPoint=Pb[TestIndex,:]
	NewTrailingPoint=CurrentPoint
	PaOrPb=2 #New current pnt is Pb

elseif Pa[TestIndex,:]!=CurrentPoint && Pa[TestIndex,:]!=TrailingPoint

	NewTriIndex=TestIndex
	NewCurrentEdge=EdgeNo
	NewCurrentPoint=Pa[TestIndex,:]
	NewTrailingPoint=CurrentPoint
	PaOrPb=1 #New current pnt is Pa

end

return NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint,PaOrPb

end


