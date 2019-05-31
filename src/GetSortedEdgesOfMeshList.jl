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
triindx=[NaN]
CurrentEdge=NaN

deleteme=[NaN NaN NaN]
FrontPointLoc=[0]
BackPointLoc=[0]
InnerPointLoc=[0]

NewTriIndex=[0]
NewCurrentEdge=[0]

NewCurrentPoint=[0. 0. 0.]
NewTrailingPoint=[0. 0. 0.]
NewInnerPoint=[0. 0. 0.]
PaOrPb=[0]

onetwo=12
twothree=23
onethree=13
NoOfChoices=[0]
rerunFunc=0
n=length(FaceNormalVector[:,1]);
removeIndx=0

#While sorted triangles still has edges we have not passed
while sum(SortedTriangles.!=0)!=length(SortedTriangles) 

	#Getting the trailing and leading pnt - finding the first triangle that is an edge in the list
	#We change edges we have passed in the list to -1 leaving only edges we have not seen in the list as '0'
	breaki=false
	for i=1:size(SortedTriangles,1) #for each tri
		for j=1:size(SortedTriangles,2) #for each potential connection

			if SortedTriangles[i,j]==0 #Its a tri with a missing connection (free edge)

				CurrentPoint=[NaN NaN NaN]# reset
				triindx=i; #First tri with free edge 
				if j==1

					CurrentEdge=12
					FirstPnt=copy(P2[triindx,:])
					TrailingPoint=copy(P1[triindx,:])
					FrontPointLoc[1]=2
					BackPointLoc[1]=1
					InnerPointLoc[1]=3

				elseif j==2

					CurrentEdge=23
					FirstPnt=copy(P3[triindx,:])
					TrailingPoint=copy(P2[triindx,:])
					FrontPointLoc[1]=3
					BackPointLoc[1]=2
					InnerPointLoc[1]=1

				elseif j==3

					CurrentEdge=13
					FirstPnt=copy(P3[triindx,:])
					TrailingPoint=copy(P1[triindx,:])
					FrontPointLoc[1]=3
					BackPointLoc[1]=1
					InnerPointLoc[1]=2

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

		#a-b is fine but going round from B is an issue as the inner triangle is not filled
		#so we can go a-b-c or a-b-d! we therefore remove tri b-c-d and rerun this
		# C• ------•                 
		#  |.＼...⁄|                  
		#  |...D•../                
		#  |. ⁄|＼•            
		# B•<  |⁄.|          
		#  |.＼•..|             
		#  |..⁄.＼|               
		# A•⁄---- •       

		if NoOfChoices[1]>1
			#We have hit a point where there are two choices
			#printstyled("Nae lookin good pal",color=:red)
			removeIndx=triindx
			rerunFunc=1 #we will need to rerun
			break
		end

		#FrontPointLoc BackPointLoc InnerPointLoc - In the current triindx which
		#points (from P1P2P3) represent the front, back and inner point
		Counter+=1; #Indx inside our lists : TrailingPoints LeadingPoints etc

		if FrontPointLoc[1]==1
			LeadingPoints[Counter,1]=triindx
		elseif FrontPointLoc[1]==2
			LeadingPoints[Counter,2]=triindx
		elseif FrontPointLoc[1]==3
			LeadingPoints[Counter,3]=triindx
		end
		if BackPointLoc[1]==1
			TrailingPoints[Counter,1]=triindx
		elseif BackPointLoc[1]==2
			TrailingPoints[Counter,2]=triindx
		elseif BackPointLoc[1]==3
			TrailingPoints[Counter,3]=triindx
		end
		if InnerPointLoc[1]==1
			InnerPoints[Counter,1]=triindx
		elseif InnerPointLoc[1]==2
			InnerPoints[Counter,2]=triindx
		elseif InnerPointLoc[1]==3
			InnerPoints[Counter,3]=triindx
		end

		if CurrentEdge==12
			SortedTriangles[triindx,1]=-1
		elseif CurrentEdge==23
			SortedTriangles[triindx,2]=-1
		elseif CurrentEdge==13
			SortedTriangles[triindx,3]=-1
		end


		#Find next outer triangle along from this triangle (they share CurrentPoint on the outer edge)
		(triindx,CurrentPoint,TrailingPoint,CurrentEdge,InnerPoint,FrontPointLoc,BackPointLoc,InnerPointLoc,NoOfChoices)=
		LoopingRoundBoundaries(triindx,CurrentPoint,TrailingPoint,CurrentEdge,
			P1,P2,P3,FeP1P2S,FeP1P3S,FeP2P3S,FrontPointLoc,BackPointLoc,InnerPointLoc,NewTriIndex,NewCurrentEdge,
								NewCurrentPoint,NewTrailingPoint,NewInnerPoint,PaOrPb,onetwo,twothree,onethree,NoOfChoices)



		if Counter>TotalNoFreeEdges
			@info Counter TotalNoFreeEdges
			error("Irk")
		end

		if FrontPointLoc[1]==BackPointLoc[1] || FrontPointLoc[1]==InnerPointLoc[1] || InnerPointLoc[1]==BackPointLoc[1]
			error("Whats happening here")
		end

	end

	if rerunFunc==1
		#Remove the bad bits
		Step=collect(1:n)
        Good=vec(fill(true,n,1))
        for i=1:length(removeIndx)
            Locs=findall(in.(Step,removeIndx[i]))
            for j=1:length(Locs)
                Good[Locs[j]]=false
            end
        end
		P1=copy(P1[Good,1:3])
        P2=copy(P2[Good,1:3])
        P3=copy(P3[Good,1:3])
        (Triangles,Points)=CreateTrianglesPointsFromP1P2P3(P1,P2,P3)
		(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
		break #exit func
	end

	if isempty(UniqueEdges)
		UniqueEdges=[1:Counter]
	else
		UniqueEdges=[UniqueEdges;[maximum(UniqueEdges[end])+1:Counter]]
	end
	CurrentPoint=[NaN NaN NaN] #reset


end

return UniqueEdges,LeadingPoints,TrailingPoints,InnerPoints,rerunFunc,P1,P2,P3,FaceNormalVector,MidPoint
end



function LoopingRoundBoundaries(triindx,CurrentPoint,TrailingPoint,CurrentEdge,
								P1,P2,P3,FeP1P2S,FeP1P3S,FeP2P3S,
								FrontPointLoc,BackPointLoc,InnerPointLoc,
								NewTriIndex,NewCurrentEdge,
								NewCurrentPoint,NewTrailingPoint,NewInnerPoint,PaOrPb,
								onetwo,twothree,onethree,NoOfChoices)
	
	#Reset values
	NewTriIndex[1]=0 #NewTriIndex
	NewCurrentEdge[1]=0 #NewCurrentEdge
	fill!(NewCurrentPoint,0.)
	fill!(NewTrailingPoint,0.)
	fill!(NewInnerPoint,0.)
	PaOrPb[1]=0
	FrontPointLoc[1]=0
	BackPointLoc[1]=0
	InnerPointLoc[1]=0

	#println("Preallocate")
	#To say we can go in two directions - if so the mesh must be cleaned before continuing
	NoOfChoices[1]=0

	InP1=ismember(P1,vec(CurrentPoint));
	if any(InP1) 
		Inside=findall(InP1)
		for i=1:length(Inside)


			if FeP1P2S.FreeFlg[Inside[i]]==true #Now check its a free edge
		
				(NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint,NewInnerPoint,FrontPointLoc,BackPointLoc,InnerPointLoc,NoOfChoices)=
				GetNewEdgePoint(Inside[i],CurrentPoint,TrailingPoint,P1,P2,P3,onetwo,
					NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint,NewInnerPoint,FrontPointLoc,BackPointLoc,InnerPointLoc,NoOfChoices)


			end
			if FeP1P3S.FreeFlg[Inside[i]]==true

				(NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint,NewInnerPoint,FrontPointLoc,BackPointLoc,InnerPointLoc,NoOfChoices)=
				GetNewEdgePoint(Inside[i],CurrentPoint,TrailingPoint,P1,P3,P2,onethree,
					NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint,NewInnerPoint,FrontPointLoc,BackPointLoc,InnerPointLoc,NoOfChoices)
				

			end

		end
	end


	InP2=ismember(P2,vec(CurrentPoint));
	if any(InP2) 
		Inside=findall(InP2)

		for i=1:length(Inside)


			if FeP1P2S.FreeFlg[Inside[i]]==true #Now check its a free edge

				(NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint,NewInnerPoint,FrontPointLoc,BackPointLoc,InnerPointLoc,NoOfChoices)=
				GetNewEdgePoint(Inside[i],CurrentPoint,TrailingPoint,P1,P2,P3,onetwo,
					NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint,NewInnerPoint,FrontPointLoc,BackPointLoc,InnerPointLoc,NoOfChoices)


			end
			if FeP2P3S.FreeFlg[Inside[i]]==true

				(NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint,NewInnerPoint,FrontPointLoc,BackPointLoc,InnerPointLoc,NoOfChoices)=
				GetNewEdgePoint(Inside[i],CurrentPoint,TrailingPoint,P2,P3,P1,twothree,
					NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint,NewInnerPoint,FrontPointLoc,BackPointLoc,InnerPointLoc,NoOfChoices)

			end

		end
	end


	InP3=ismember(P3,vec(CurrentPoint));
	if any(InP3) 
		Inside=findall(InP3)
		for i=1:length(Inside)

			if FeP1P3S.FreeFlg[Inside[i]]==true #Now check its a free edge

				(NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint,NewInnerPoint,FrontPointLoc,BackPointLoc,InnerPointLoc,NoOfChoices)=
				GetNewEdgePoint(Inside[i],CurrentPoint,TrailingPoint,P1,P3,P2,onethree,
					NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint,NewInnerPoint,FrontPointLoc,BackPointLoc,InnerPointLoc,NoOfChoices)

			end
			if FeP2P3S.FreeFlg[Inside[i]]==true

				(NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint,NewInnerPoint,FrontPointLoc,BackPointLoc,InnerPointLoc,NoOfChoices)=
				GetNewEdgePoint(Inside[i],CurrentPoint,TrailingPoint,P2,P3,P1,twothree,
					NewTriIndex,NewCurrentEdge,NewCurrentPoint,NewTrailingPoint,NewInnerPoint,FrontPointLoc,BackPointLoc,InnerPointLoc,NoOfChoices)


			end

		end
	end	

	triindx=NewTriIndex;
	CurrentPoint=NewCurrentPoint;
	TrailingPoint=NewTrailingPoint;
	CurrentEdge=NewCurrentEdge;
	InnerPoint=NewInnerPoint


return triindx,CurrentPoint,TrailingPoint,CurrentEdge,InnerPoint,FrontPointLoc,BackPointLoc,InnerPointLoc,NoOfChoices

end


function GetNewEdgePoint(TestIndex,CurrentPoint,TrailingPoint,
						Pa,Pb,Pc,EdgeNo,NewTriIndex,
						NewCurrentEdge,NewCurrentPoint,NewTrailingPoint,NewInnerPoint,
						FrontPointLoc,BackPointLoc,InnerPointLoc,NoOfChoices)

if Pb[TestIndex,:]!=CurrentPoint && Pb[TestIndex,:]!=TrailingPoint

	NewTriIndex=TestIndex
	NewCurrentEdge=EdgeNo
	NewCurrentPoint=Pb[TestIndex,:]
	NewTrailingPoint=CurrentPoint
	NewInnerPoint=copy(Pc[TestIndex,:])

	if EdgeNo==12
		FrontPointLoc[1]=2;
		BackPointLoc[1]=1;
		InnerPointLoc[1]=3;
	elseif EdgeNo==13
		FrontPointLoc[1]=3;
		BackPointLoc[1]=1;
		InnerPointLoc[1]=2;
	elseif EdgeNo==23
		FrontPointLoc[1]=3;
		BackPointLoc[1]=2;
		InnerPointLoc[1]=1;		
	end
	NoOfChoices[1]+=1

elseif Pa[TestIndex,:]!=CurrentPoint && Pa[TestIndex,:]!=TrailingPoint

	NewTriIndex=TestIndex
	NewCurrentEdge=EdgeNo
	NewCurrentPoint=Pa[TestIndex,:]
	NewTrailingPoint=CurrentPoint
	NewInnerPoint=copy(Pc[TestIndex,:])

	if EdgeNo==12		
		FrontPointLoc[1]=1;
		BackPointLoc[1]=2;
		InnerPointLoc[1]=3;
	elseif EdgeNo==13
		FrontPointLoc[1]=1;
		BackPointLoc[1]=3;
		InnerPointLoc[1]=2;
	elseif EdgeNo==23
		FrontPointLoc[1]=2;
		BackPointLoc[1]=3;
		InnerPointLoc[1]=1;		
	end
	NoOfChoices[1]+=1

end

return NewTriIndex,NewCurrentEdge,
NewCurrentPoint,NewTrailingPoint,NewInnerPoint,
FrontPointLoc,BackPointLoc,InnerPointLoc,NoOfChoices

end


