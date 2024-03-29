function CheckIfInsidePreviousBoundary(MidPoint,P1,P2,P3,FaceNormalVector,
                P1New,P2New,P3New,MidPointNew,max_target_edge_length)

#Previous mesh:
#MidPoint,P1,P2,P3,FaceNormalVector

#New mesh:
#P1New,P2New,P3New

#=
using Plots
fig = plot()
scatter!([MidPointNew[:,1]],[MidPointNew[:,3]])
scatter!([MidPoint[:,1]],[MidPoint[:,3]])
=#


#Get edge triangles
(P1P2FreeFlg,P2P3FreeFlg,P1P3FreeFlg)=EdgeConstraints(P1New,P2New,P3New,MidPointNew);
#Edge points of the new mesh [x y z]
p1p2idx=findall(P1P2FreeFlg)
p1p3idx=findall(P1P3FreeFlg)
p2p3idx=findall(P2P3FreeFlg)
EdgePnts=[P1New[p1p2idx,:];
          P2New[p1p2idx,:];
          P1New[p1p3idx,:];
          P3New[p1p3idx,:];
          P2New[p2p3idx,:];
          P3New[p2p3idx,:]];

#Find boundary of the previous mesh:
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

#New coords
V1=[0. 0. 1. ]; #Pointing up
#List we check for each edge if contained within our set padding
PntList=fill(false,length(EdgePnts[:,1]))

#n=length(EdgePnts[:,1])
#X=zeros(n);
#Y=zeros(n);
#Z=zeros(n);

#Radius based on a 1/4 max bound
(Tris,Pnts)=CutAndDisplaceJulia.CreateTrianglesPointsFromP1P2P3(P1,P2,P3)
maxX=maximum(Pnts[:,2])
minX=minimum(Pnts[:,2])
maxY=maximum(Pnts[:,3])
minY=minimum(Pnts[:,3])
maxZ=maximum(Pnts[:,4])
minZ=minimum(Pnts[:,4])
Ranges=maximum([maxX-minX maxY-minY maxZ-minZ]);
SetRadius=Ranges/10 #10% of max range (any dir). max_target_edge_length

lps=0
#For each edge loop - here we check if any points on the new boundary are outside a set cylinder radius of the old boundary
for i=1:length(UniqueEdges)

  b=vec(UniqueEdges[i])

  for j=b

    #Extract the points on the current bit of the edge
    (Pa,~) =GrabPointNew3(LeadingPoints,P1,P2,P3,j)  #LeadingPoint
    (Pb,~) =GrabPointNew3(TrailingPoints,P1,P2,P3,j) #TrailingPoint
    #Midpoint of edge
    Pc=  [((Pa[1]+Pb[1])./2) ((Pa[2]+Pb[2])./2) ((Pa[3]+Pb[3])./2)];
    #Vector along the edge of the old mesh
    v=normr([(Pa[1]-Pb[1]) (Pa[2]-Pb[2]) (Pa[3]-Pb[3])]);
    #Length between Pa and Pb
    l=sqrt(((Pa[1]-Pb[1])^2)+((Pa[2]-Pb[2])^2)+((Pa[3]-Pb[3])^2));
    
    X=copy(EdgePnts[:,1]);
    Y=copy(EdgePnts[:,2]);
    Z=copy(EdgePnts[:,3]);

    #Rotate all new edge points so the centre is Pc and PaPb is the Z-axis
    (X,Y,Z) =RotateObject3DAllignVectors(v,V1,X,Y,Z,Pc[1],Pc[2],Pc[3]);

    #put midpoint at 0,0,0
    X=X.-Pc[1];
    Y=Y.-Pc[2];
    Z=Z.-Pc[3];

    (~,ρ,z)=cart2pol(X,Y,Z)

    #If all
    zbounds=(l/2+(max_target_edge_length/2)); #half the total length with some extra padding at either end
    for j=1:length(ρ)
      if ρ[j]<SetRadius && abs(z[j])<zbounds
        PntList[j]=true
      end
    end

    end

end

@bp 

if all(PntList)
  InsidePrevious=true
else
  InsidePrevious=false
end

return(InsidePrevious)

end


function GrabPointNew3(PointsIdxList,P1,P2,P3,j)
#Extract the points on the current bit of the edge
Point=[0. 0. 0.]
Indx=0
for k=1:3
    Idx=PointsIdxList[j,k]
    if Idx==0
        continue
    elseif k==1
        Point=P1[Idx,:]
        Indx=Idx
    elseif k==2
        Point=P2[Idx,:]
        Indx=Idx
    elseif k==3
        Point=P3[Idx,:]
        Indx=Idx
    end
end

return Point,Indx
end

