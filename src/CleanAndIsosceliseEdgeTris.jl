function CleanAndIsosceliseEdgeTris(MidPoint,P1,P2,P3,Triangles,FaceNormalVector)

    ## Firstly we collapse edge triangles that share an inner point
    #This loop assumes these double edge bits only have two connections
    #
    #  ¯.¯\¯¯|¯¯/¯.¯       ¯.¯\¯¯ ¯¯/¯.¯
    #  ....\ | /....  -- > ....\   /....
    #  .....\|/.....       .....\ /.....
    #  .............       .............      
(newTris,removeIndx)=CutAndDisplaceJulia.CollapseEdgeTris(P1,P2,P3,MidPoint,FaceNormalVector)
n=length(Triangles[:,1]);




## Get edge triangles
println("Grabbing edge tris") 
#Get edge triangles
(P1P2FreeFlg,P2P3FreeFlg,P1P3FreeFlg)=EdgeConstraints(P1,P2,P3,MidPoint);
(SortedTriangles,ConnectedEdge)=ConnectedConstraints(P1,P2,P3,MidPoint)

#n*3 list of triangles connected to this tri





#Only doing if there are changes
if removeIndx!=[0 0] || newTris!=[0 0 0 0 0 0 0 0 0]

    #Remove first dud value
    newTris=copy(newTris[2:end,:])
    removeIndx=copy(removeIndx[2:end])

    
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
    Points=zeros(Int(length(P1)),3)
    Points[1:3:end,:]=copy(P1);
    Points[2:3:end,:]=copy(P2);
    Points[3:3:end,:]=copy(P3);
    n=length(Points)
    Points=[1:n/3 Points]
    Triangles=fill(0,Int(n/9),3)
    Triangles[:,1]=1:3:n/3;
    Triangles[:,2]=2:3:n/3;
    Triangles[:,3]=3:3:n/3;


    #(P1,P2,P3) = CreateP1P2P3( Triangles,Points ); 
    try (FaceNormalVector,MidPoint) = CreateFaceNormalAndMidPoint(Points,Triangles)
    catch 
        println("Check your surface, more than 2 duplicate edge tris?")
        error("Remesh here")
    end
    ### Get edge triangles (Of cleaned tri)
    #(P1P2FreeFlg,P2P3FreeFlg,P1P3FreeFlg)=EdgeConstraints(P1,P2,P3,MidPoint);

end

## PART 2: Now Rotate so always isosceles tris on edge
(FeP1P2S,FeP1P3S,FeP2P3S)=GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector)

#Do for P1 P2: (Function at base of file)
(P1,P2,P3)=MakeEqEdgeTris(FeP1P2S,P1,P2,P3,MidPoint,FaceNormalVector); #

#Do for P1 P3: (Function at base of file)
(P1,P2,P3)=MakeEqEdgeTris(FeP1P3S,P1,P3,P2,MidPoint,FaceNormalVector);

#Do for P2 P3: (Function at base of file)
(P1,P2,P3)=MakeEqEdgeTris(FeP2P3S,P2,P3,P1,MidPoint,FaceNormalVector);


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

#=
FeM2Ev=Array{Float64,2}(undef, n,3);FeM2Ev=fill!(FeM2Ev, NaN)
FeM2Ev[Indx,:]=normr([T.FeMd[Indx,1]-MidPoint[Indx,1] T.FeMd[Indx,2]-MidPoint[Indx,2] T.FeMd[Indx,3]-MidPoint[Indx,3]]);
FeEv=Array{Float64,2}(undef, n,3);FeEv=fill!(FeEv, NaN)
FeEv[Indx,:]=normr([(Pa[Indx,1]-Pb[Indx,1]) (Pa[Indx,2]-Pb[Indx,2]) (Pa[Indx,3]-Pb[Indx,3])]);
##
#Fix to make sure that the the Edge vector is counter clockwise from the
#mid2edge vector when looking in the normal direction:
for i=1:length(Indx)
    #First rotate to flat:
    V1=[0 0 1]; #Pointing up
    #Get Nx Ny Nz for the vectors and put in here. 
    X=[FeM2Ev[Indx[i],1] FeEv[Indx[i],1]];
    Y=[FeM2Ev[Indx[i],2] FeEv[Indx[i],2]];
    Z=[FeM2Ev[Indx[i],3] FeEv[Indx[i],3]];

    #Rotate so vectors are flat:
    (X,Y,~) = RotateObject3DAllignVectors(FaceNormalVector[Indx[i],:],V1,X,Y,Z,0,0,0);
    #Now get the two vectors as 2D coords (2nd we rotate by 90 counter
    #clockwise)
    VM2Ev=[X[1] Y[1]]; 
    VEv=[Y[2] -X[2]]; 
    #Calculate the dot product
    AllignFlag=dot(VM2Ev,VEv);
    #See if these allign or not:
    if AllignFlag>0
        FeEv[Indx[i],:]=-FeEv[Indx[i],:];
    end
end
=#

#Internal angles
Ang=fill(NaN,n,1)
#First we recompute the mid2edvec length (perp):
#Angle between vectors: 

@info size(T.FeEv) size(FeIn2Ev) Indx

Upsidedown=FaceNormalVector[Indx,3].<0;
for i=1:length(Indx)
    idx=Indx[i];
    Ang[idx]=(pi/2)-(acos(dot(vec(FeIn2Ev[idx,:]),vec(T.FeEv[idx,:]))));
    if Upsidedown[i]
        Ang[idx]=-Ang[idx];
    end
end
# #Length of R
# FePc2ELe(Flag,:);
# #Centre of Rotation
# FeMd(Flag,:);
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

