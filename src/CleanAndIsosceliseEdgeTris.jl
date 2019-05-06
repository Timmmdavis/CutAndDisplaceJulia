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

error("heere")
##
println("Collate") 

Points=[zeros(size(P1));zeros(size(P1));zeros(size(P1))];
Points[1:3:end]=P1;
Points[2:3:end]=P2;
Points[3:3:end]=P3;
Triangles=1:1:length(Points)/3;
Triangles=reshape(Triangles,3,[])';

#figure;
#trisurf(Triangles,Points(:,1),Points(:,2),Points(:,3));
#Points=[(1:1:length(Points(:,1)))",Points];

## PART 2: Now Rotate so always eq lat tris on edge
#[P1,P2,P3] = CreateP1P2P3( Triangles,Points ); 
(MidPoint,FaceNormalVector) = MidPointCreate(Points,Triangles,0);

## Get edge triangles (Of cleaned tri)
(P1P2FreeFlg,P2P3FreeFlg,P1P3FreeFlg)=EdgeCons(P1,P2,P3,MidPoint);

#hold on
#scatter3(MidPoint(P1P2FreeFlg,1),MidPoint(P1P2FreeFlg,2),MidPoint(P1P2FreeFlg,3),"filled","red")
#scatter3(MidPoint(P2P3FreeFlg,1),MidPoint(P2P3FreeFlg,2),MidPoint(P2P3FreeFlg,3),"filled","green")
#scatter3(MidPoint(P1P3FreeFlg,1),MidPoint(P1P3FreeFlg,2),MidPoint(P1P3FreeFlg,3),"filled","blue")
#title("CleanedDupEdges - FreeEdgeIndx - Going into Equi eqs")

## Get new point if all tris are now isos
#Do for P1 P2: (Function at base of file)
(P1,P2,P3)=MakeEqEdgeTris(P1P2FreeFlg,P1,P2,P3,MidPoint,FaceNormalVector);

#Do for P1 P3: (Function at base of file)
(P1,P2,P3)=MakeEqEdgeTris(P1P3FreeFlg,P1,P3,P2,MidPoint,FaceNormalVector);

#Do for P2 P3: (Function at base of file)
(P1,P2,P3)=MakeEqEdgeTris(P2P3FreeFlg,P2,P3,P1,MidPoint,FaceNormalVector);


Points=[zeros(size(P1));zeros(size(P1));zeros(size(P1))];
Points[1:3:end]=P1;
Points[2:3:end]=P2;
Points[3:3:end]=P3;
Triangles=1:1:length(Points)/3;
Triangles=reshape(Triangles,3,[])';

#figure;
#trisurf(Triangles,Points(:,1),Points(:,2),Points(:,3));
#axis("equal")

Points=[(1:1:length(Points(:,1)))',Points];
(MidPoint,FaceNormalVector) = MidPointCreate(Points,Triangles,0);

return P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector

end




function MakeEqEdgeTris(Flg,Pa,Pb,Pc,MidPoint,FaceNormalVector)
#Fills array with values if the connection is a free edge. 

#Initialise some arrays:

#Lengths of Free edges
FeLe=nan(length(MidPoint)/3,1);
#Length from innerpoint to edge Midpoint
FeIn2ELe=nan(length(MidPoint)/3,1);
#MidPoints of Free edges
FeMd=nan(length(MidPoint)/3,3);
#Vector parallel to Free edges (edge vector Ev)
FeEv=nan(length(MidPoint)/3,3);
#Vector pointing from midpoint to midpoint of Free edge (mid to edge vector
#M2Ev)
FeM2Ev=nan(length(MidPoint)/3,3);
#Internal angles
IntAng=nan(length(MidPoint)/3,1);
#Vector from inner point to edge midpoint
FeIn2Ev=nan(length(MidPoint)/3,3);

#Create index thats says if the connection is a free edge. 
Flag=logical(Flg); 
#Midpoint of edge
FeMd[Flag,:]=  [((Pa[Flag,1]+Pb[Flag,1])/2), ((Pa[Flag,2]+Pb[Flag,2])/2), ((Pa[Flag,3]+Pb[Flag,3])/2)];
#Vector from midpoint to edge midpoint. 
FeM2Ev[Flag,:]=normr([FeMd[Flag,1]-MidPoint[Flag,1],FeMd[Flag,2]-MidPoint[Flag,2],FeMd[Flag,3]-MidPoint[Flag,3]]);
#Vector from inner point to edge midpoint
FeIn2Ev[Flag,:]=normr([FeMd[Flag,1]-Pc[Flag,1],FeMd[Flag,2]-Pc[Flag,2],FeMd[Flag,3]-Pc[Flag,3]]);

#Vector pointing along edge
FeEv[Flag,:]=normr([(Pa[Flag,1]-Pb[Flag,1]),  (Pa[Flag,2]-Pb[Flag,2]),     (Pa[Flag,3]-Pb[Flag,3])]);

#Internal angle of the triangle (angle between edges that are not the free
#edge in question). 
v=normr([(Pa[Flag,1]-Pc[Flag,1]), (Pa[Flag,2]-Pc[Flag,2]), (Pa[Flag,3]-Pc[Flag,3])]);
w=normr([(Pb[Flag,1]-Pc[Flag,1]), (Pb[Flag,2]-Pc[Flag,2]), (Pb[Flag,3]-Pc[Flag,3])]);

Indx=find(Flag);
for i = 1:length(v(:,1))
    IntAng[Indx[i]]=rad2deg(acos(dot(v[i,:],w[i,:])));
end
#|a| is the magnitude (length) of vector a
#|b| is the magnitude (length) of vector b
#? is the angle between a and b
#a · b = |a| × |b| × cos(?) 
#so acos the dot of the normalised gives theta


##
#Fix to make sure that the the Edge vector is counter clockwise from the
#mid2edge vector when looking in the normal direction:
for i=1:length(Indx)
    #First rotate to flat:
    V1=[0,0,1]; #Pointing up
    #Get Nx Ny Nz for the vectors and put in here. 
    X=[FeM2Ev[Indx[i],1],FeEv[Indx[i],1]];
    Y=[FeM2Ev[Indx[i],2],FeEv[Indx[i],2]];
    Z=[FeM2Ev[Indx[i],3],FeEv[Indx[i],3]];
    #Rotate so vectors are flat:
    (X,Y,~) = RotateObject3dAllignVectors(FaceNormalVector[Indx[i],:],V1,X,Y,Z,0,0,0);
    #Now get the two vectors as 2D coords (2nd we rotate by 90 counter
    #clockwise)
    VM2Ev=[X[1], Y[1]]; 
    VEv=[Y[2], -X[2]]; 
    #Calculate the dot product
    AllignFlag=dot(VM2Ev,VEv);
    #See if these allign or not:
    if AllignFlag>0
        FeEv[Indx[i],:]=-FeEv[Indx[i],:];
    end
end

## 
#If triangles are not completly equilateral then the mid2Ed vec calculated
#above is not perpendicular to the edge (meaning poor calculation of K2).
#This fixes this for non eq tris.


#First we recompute the mid2edvec length (perp):
#Angle between vectors: 
Ang=pi/2-(acos(dot(FeIn2Ev[Flag,:]',FeEv[Flag,:]')));
Upsidedown=(FaceNormalVector[Flag,3]<0)';
Ang[Upsidedown]=-Ang[Upsidedown];
# #Length of R
# FePc2ELe(Flag,:);
# #Centre of Rotation
# FeMd(Flag,:);
#Default axis
Vect=FaceNormalVector[Flag,:];
#Place Point to be rotated in correct pos
PcCoords=Pc[Flag,:]-FeMd[Flag,:];
# EXAMPLE:
#     Rotate point (1;2;3) around vector (4;5;6) by an angle of pi/2
#     P = [1;2;3];  # create the point
#     V = [4;5;6];  # create vector around which rotation is performed
#     Qrot = qGetRotQuaternion( pi/2, V );
#     P2 = qRotatePoint( P, Qrotate ); ];
#
Pc2=zeros(length(Ang),3);
for i=1:length(Ang)
    Qrot = qGetRotQuaternion( Ang[i], Vect[i,:]' );
    Pc2[i,:] = qRotatePoint( PcCoords[i,:]', Qrot ); 
end

#Move back to orig coords:
Pc2=Pc2+FeMd[Flag,:];


# ## Only if you want eq tris
# #Length of edge
# FeLe(Flag,:)=sqrt(((Pa(Flag,1)-Pb(Flag,1)).^2)+((Pa(Flag,2)-Pb(Flag,2)).^2)+((Pa(Flag,3)-Pb(Flag,3)).^2));
# #Length of innerpoint to edge innerpoint
# FeIn2ELe(Flag,:)=sqrt(((FeMd(Flag,1)-Pc(Flag,1)).^2)+((FeMd(Flag,2)-Pc(Flag,2)).^2)+((FeMd(Flag,3)-Pc(Flag,3)).^2));
# Shdbe=(FeLe*sqrt(3))/2;
# ToMove=Shdbe-FeIn2ELe;
# FeIn2Ev(Flag,:)=normr([FeMd(Flag,1)-Pc2(:,1),FeMd(Flag,2)-Pc2(:,2),FeMd(Flag,3)-Pc2(:,3)]);
# Pc2=Pc2+(-ToMove(Flag,:).*FeIn2Ev(Flag,:));
# ##


#hold on
#scatter3(Pc2(:,1),Pc2(:,2),Pc2(:,3),"g","filled")

#And clean up connected tris
for i=1:length(Indx)
    InPa=ismember(Pa(:,:),Pc(Indx(i),:),"rows");
    InPb=ismember(Pb(:,:),Pc(Indx(i),:),"rows");
    InPc=ismember(Pc(:,:),Pc(Indx(i),:),"rows");
    InPaIndx=find(InPa);
    for j=1:length(InPaIndx)
        Pa[InPaIndx(j),:]=Pc2[i,:];
    end
    
    InPbIndx=find(InPb);
    for j=1:length(InPbIndx)      
        Pb[InPbIndx(j),:]=Pc2[i,:];
    end    
    
    InPcIndx=find(InPc);
    for j=1:length(InPcIndx)    
        Pc[InPcIndx(j),:]=Pc2[i,:];
    end    
end

Pc[Flag,:]=Pc2;

return Pa,Pb,Pc

end

