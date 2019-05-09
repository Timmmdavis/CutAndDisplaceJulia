function CleanAndIsosceliseEdgeTris(MidPoint,P1,P2,P3,Triangles,FaceNormalVector)

## Get edge triangles
println("Grabbing edge tris") 
#Get edge triangles
(P1P2FreeFlg,P2P3FreeFlg,P1P3FreeFlg)=EdgeConstraints(P1,P2,P3,MidPoint);
(SortedTriangles,ConnectedEdge)=ConnectedConstraints(P1,P2,P3,MidPoint)

#n*3 list of triangles connected to this tri


## Firstly we collapse edge triangles that share an inner point
#This loop assumes these double edge bits only have two connections
#
#  ¯.¯\¯¯|¯¯/¯.¯       ¯.¯\¯¯ ¯¯/¯.¯
#  ....\ | /....  -- > ....\   /....
#  .....\|/.....       .....\ /.....
#  .............       .............      

n=length(Triangles[:,1]);
(SixPntsP1P2)=CreateSortedEdgePoints(P1,P2);
(SixPntsP2P3)=CreateSortedEdgePoints(P2,P3);
(SixPntsP3P1)=CreateSortedEdgePoints(P3,P1);
#number of connected tris (Sorted tris rows not == to 0)
NoConnections=sum(SortedTriangles.!=0,dims=2)
EdgeTri=NoConnections.<3 #- edge tri

fillme=zeros(1,9)
fillme2=fill(0,1,2)
newTris=zeros(1,9)
removeIndx=fill(0,1,2);
AvgVect=zeros(1,3);
append=0;
for i=1:n #for each triangle
    
    if EdgeTri[i] #current triangle is an edge
        #We look at each connection of this tri
        for j=1:3
            Next2EdgeTri=SortedTriangles[i,j]
            if Next2EdgeTri==0
                continue
            end

            if EdgeTri[Next2EdgeTri] #both edges and connected
                
                #No point doing this for the 2nd tri
                if any(in.(i,removeIndx))
                    continue
                end
                

                #Getting connected edges of each tri (P1 P2 or P3)
                #j defines if its Col1=P1P2 Col2=P2P3 Col3=P1P3
                if j==1
                    #EdgeTriEdge=12
                    EdgeTriNonConnectedPoint=P3[i,:]
                    if P2P3FreeFlg[i]
                        EdgeTriInnerPoint=P1[i,:]
                    else
                        EdgeTriInnerPoint=P2[i,:]
                    end
                elseif j==2
                    #EdgeTriEdge=23
                    EdgeTriNonConnectedPoint=P1[i,:]
                    if P1P2FreeFlg[i]
                        EdgeTriInnerPoint=P3[i,:]
                    else
                        EdgeTriInnerPoint=P2[i,:]
                    end                    
                elseif j==3
                    #EdgeTriEdge=13
                    EdgeTriNonConnectedPoint=P2[i,:]
                    if P1P2FreeFlg[i]
                        EdgeTriInnerPoint=P3[i,:]
                    else
                        EdgeTriInnerPoint=P1[i,:]
                    end                       
                end

                Next2EdgeTriEdge=ConnectedEdge[i,j]
                if Next2EdgeTriEdge==12
                    Next2EdgeTriNonConnectedPoint=P3[Next2EdgeTri,:]
                    #if P2P3FreeFlg[Next2EdgeTri]
                    #    Next2EdgeTriInnerPoint=P1[Next2EdgeTri,:]
                    #else
                    #    Next2EdgeTriInnerPoint=P2[Next2EdgeTri,:]
                    #end                    
                elseif Next2EdgeTriEdge==23
                    Next2EdgeTriNonConnectedPoint=P1[Next2EdgeTri,:]
                    #if P1P2FreeFlg[Next2EdgeTri]
                    #    Next2EdgeTriInnerPoint=P3[Next2EdgeTri,:]
                    #else
                    #    Next2EdgeTriInnerPoint=P2[Next2EdgeTri,:]
                    #end          
                elseif Next2EdgeTriEdge==13
                    Next2EdgeTriNonConnectedPoint=P2[Next2EdgeTri,:]
                    #if P1P2FreeFlg[Next2EdgeTri]
                    #    Next2EdgeTriInnerPoint=P3[Next2EdgeTri,:]
                    #else
                    #    Next2EdgeTriInnerPoint=P1[Next2EdgeTri,:]
                    #end   
                end

                

                #should match - Next2EdgeTriInnerPoint EdgeTriInnerPoint
               
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
        end
    end
end

#Only doing if there are changes
if removeIndx!=[0 0] || newTris!=[0 0 0 0 0 0 0 0 0]

    #Now redefine P1 P2 P3 
    n_new=n-length(removeIndx)+size(newTris,1)
    P1new=zeros(n_new,3)
    P2new=zeros(n_new,3)
    P3new=zeros(n_new,3)

    #should always be more than new tris
    counter=1;
    for i=1:n
        if any(in.(i,removeIndx))
            #skip
        else
            P1new[counter,:]=P1[i,:]
            P2new[counter,:]=P2[i,:]
            P3new[counter,:]=P3[i,:]
            counter+=1
        end
    end
    P1new=[P1new;newTris[:,1:3]]
    P2new=[P2new;newTris[:,4:6]]
    P3new=[P3new;newTris[:,7:9]]


    P1=P1new
    P2=P2new
    P3=P3new

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
    #(P1,P2,P3) = CreateP1P2P3( Triangles,Points ); 
    try (MidPoint,FaceNormalVector) = CreateFaceNormalAndMidPoint(Points,Triangles)
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

