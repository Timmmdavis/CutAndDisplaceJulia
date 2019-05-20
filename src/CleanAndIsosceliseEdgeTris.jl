function CleanAndIsosceliseEdgeTris(MidPoint,P1,P2,P3,Triangles,FaceNormalVector)

#Remove slither tris
(SortedTriangles,ConnectedEdge)=ConnectedConstraints(P1,P2,P3,MidPoint)
#number of connected tris (Sorted tris rows not == to 0)
NoConnections=sum(SortedTriangles.!=0,dims=2)

#Remove any slither tris
( Area,HalfPerimeter ) = CutAndDisplaceJulia.AreaOfTriangle3D( P1[:,1],P1[:,2],P1[:,3],P2[:,1],P2[:,2],P2[:,3],P3[:,1],P3[:,2],P3[:,3] );
Good=vec(fill(true,length(Area)))
tol=mean(Area)/3
for i=1:length(Good)
    if Area[i]<tol && NoConnections[i]<2 
        Good[i]=false
    end
end

SlithersRemoved=sum(Good.==false)
println("$SlithersRemoved slither triangles have been removed")
P1=copy(P1[Good,1:3])
P2=copy(P2[Good,1:3])
P3=copy(P3[Good,1:3])
(Points,Triangles)=CreateTrianglesPointsFromP1P2P3(P1,P2,P3)
(FaceNormalVector,MidPoint) = CreateFaceNormalAndMidPoint(Points,Triangles)

rerunFunc=1 #sometime we need to run twice
while rerunFunc==1

    (newTris,removeIndx,rerunFunc)=CutAndDisplaceJulia.CollapseEdgeTris(P1,P2,P3,MidPoint,FaceNormalVector)
    n=length(Triangles[:,1]);


    #Get edge triangles
    (P1P2FreeFlg,P2P3FreeFlg,P1P3FreeFlg)=EdgeConstraints(P1,P2,P3,MidPoint);
    (SortedTriangles,ConnectedEdge)=ConnectedConstraints(P1,P2,P3,MidPoint);

    #n*3 list of triangles connected to this tri


    #Only doing if there are changes
    if size(removeIndx)!=()

        #Remove first dud value
        newTris=copy(newTris[2:end,:])
        removeIndx=copy(removeIndx[2:end])

        NoCollapsed=length(removeIndx)
        NoCreated=size(newTris,1)
        if NoCollapsed>0
            println("$NoCollapsed shared inner point edge triangles have been collapsed into $NoCreated triangles with unique inner points")
        end

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

        P1=[P1;newTris[:,1:3]]
        P2=[P2;newTris[:,4:6]]
        P3=[P3;newTris[:,7:9]]


        ## Recreate tri
        (Points,Triangles)=CreateTrianglesPointsFromP1P2P3(P1,P2,P3)
        try (FaceNormalVector,MidPoint) = CreateFaceNormalAndMidPoint(Points,Triangles)
        catch 
            println("Check your surface, more than 2 duplicate edge tris?")
            error("Remesh here")
        end
        ### Get edge triangles (Of cleaned tri)
        #(P1P2FreeFlg,P2P3FreeFlg,P1P3FreeFlg)=EdgeConstraints(P1,P2,P3,MidPoint);

    end
end

# The aim is to isoscelise edge triangles by swinging the inner point around
# the edge midpoint
#
#                       Em
# ¯.¯\¯¯¯¯¯ ⁄¯.¯  ¯.¯\¯¯•¯¯ ⁄¯.¯      ¯.¯\¯¯ ¯¯/¯.¯
# ....\   ⁄.....  ....\   ⁄..... -- > ....\   /....
# .....\⁄.......  .....\⁄.......      .....\ /.....
# .....I........  .|...I........       .....I.......     
#                  \  
#                   ＼ _ >
#                          
#  
# Em = Edge MidPoint
# I=Inner Point

## PART 2: Now Rotate so always isosceles tris on edge
(FeP1P2S,FeP1P3S,FeP2P3S)=GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector);

#Do for P1 P2: (Function at base of file)
(P1,P2,P3)=MakeEqEdgeTris(FeP1P2S,P1,P2,P3,MidPoint,FaceNormalVector); #

#Do for P1 P3: (Function at base of file)
(P1,P3,P2)=MakeEqEdgeTris(FeP1P3S,P1,P3,P2,MidPoint,FaceNormalVector);

#Do for P2 P3: (Function at base of file)
(P2,P3,P1)=MakeEqEdgeTris(FeP2P3S,P2,P3,P1,MidPoint,FaceNormalVector);


## Recreate tri
Points=zeros(Int(length(P1)),3)
Points[1:3:end,:]=P1;
Points[2:3:end,:]=P2;
Points[3:3:end,:]=P3;
n=length(Points)
Points=[1:n/3 Points]
Triangles=fill(0,Int(n/9),3)
Triangles[:,1]=1:3:n/3;
Triangles[:,2]=2:3:n/3;
Triangles[:,3]=3:3:n/3;

#[P1,P2,P3] = CreateP1P2P3( Triangles,Points ); 
(FaceNormalVector,MidPoint) = CreateFaceNormalAndMidPoint(Points,Triangles)


return P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector

end




function MakeEqEdgeTris(T,Pa,Pb,Pc,MidPoint,FaceNormalVector)
#Fills array with values if the connection is a free edge. 

#= Structure contains
TriangleEdges (T) =
FeLe;FeMd;FeEv;FeM2Ev;FreeFlg;FeM2ELe;IntAng;
all defined in: GetCrackTipElements3D
=#

## 
#If triangles are not completly equilateral then the mid2Ed vec calculated
#above is not perpendicular to the edge (meaning poor calculation of K2).
#This fixes this for non eq tris.

n=Int(length(MidPoint)/3);
Indx=findall(T.FreeFlg);

#Vector from inner point to edge midpoint

FeIn2Ev=Array{Float64,2}(undef, n,3);FeIn2Ev=fill!(FeIn2Ev, NaN)
FeIn2Ev[Indx,:]=normr([T.FeMd[Indx,1]-Pc[Indx,1] T.FeMd[Indx,2]-Pc[Indx,2] T.FeMd[Indx,3]-Pc[Indx,3]]);

#Internal angles
Ang=fill(NaN,n,1)
#First we recompute the mid2edvec length (perp):
#Angle between vectors: 


#Upsidedown=FaceNormalVector[Indx,3].>0;
for i=1:length(Indx)
    idx=Indx[i];
    Ang[idx]=(pi/2)-(acos(dot(vec(FeIn2Ev[idx,:]),vec(T.FeEv[idx,:]))));

end

#Default axis
Vect=FaceNormalVector[Indx,:];
#Place Point to be rotated in correct pos
PcCoords=Pc[Indx,:].-T.FeMd[Indx,:];
# EXAMPLE:
#     Rotate point (1;2;3) around vector (4;5;6) by an angle of pi/2
#     P = [1;2;3];  # create the point
#     V = [4;5;6];  # create vector around which rotation is performed
#     Qrot = qGetRotQuaternion( pi/2, V );
#     P2 = qRotatePoint( P, Qrotate ); ];
#
Pc2=zeros(length(Indx),3);
for i=1:length(Indx)
    idx=Indx[i];    
    Qrot = qGetRotQuaternion(  Ang[idx], vec(Vect[i,:]) );
    Pc2[i,:] = qRotatePoint( PcCoords[i,:], Qrot ); 
end

#Move back to orig coords:
Pc2=Pc2.+T.FeMd[Indx,:];


#And clean up connected tris
for i=1:length(Indx)

    OldPoint=Pc[Indx[i],:]
    NewPoint=Pc2[i,:];

    #if the connected tri contains the old point
    InPa=ismember(Pa,OldPoint);
    InPb=ismember(Pb,OldPoint);
    InPc=ismember(Pc,OldPoint);

    #loop through and update to new point
    for j=1:length(InPa)
        if InPa[j]
            Pa[j,:]=NewPoint;
        end
    end
    
    for j=1:length(InPb)
        if InPb[j]
            Pb[j,:]=NewPoint;
        end
    end    
    
    for j=1:length(InPc)  
        if InPc[j]
            Pc[j,:]=NewPoint;
        end
    end    
end

Pc[Indx,:]=Pc2;

return Pa,Pb,Pc

end

