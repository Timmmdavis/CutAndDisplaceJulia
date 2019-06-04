function IsosceliseEdgeTrisNew(MidPoint,P1,P2,P3,Triangles,Points,FaceNormalVector)

#=
##################################################################
println("Remove me and add 2 inputs:")
(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
##################################################################
=#

(UniqueEdges,LeadingPoints,TrailingPoints,InnerPoints,rerunFunc,P1,P2,P3,FaceNormalVector,MidPoint)=
CutAndDisplaceJulia.GetSortedEdgesOfMeshList(P1,P2,P3,FaceNormalVector,MidPoint)

(SortedTriangles,ConnectedEdge)=CutAndDisplaceJulia.ConnectedConstraints(P1,P2,P3,MidPoint);

# The aim is to isoscelise edge triangles by finding the intersection between the edge MidPointVec and 
# the longest inner edge. We then retriangulate with x as the new inner point
#
#   q                   Em             q
# ¯.¯\¯¯¯¯¯⁄¯.¯   ¯.¯\¯¯¯•¯⁄¯.¯      ¯.¯\¯•¯⁄¯.¯
# ....\   ⁄.....  ....\  ↓⁄..... -- > ...\ ⁄....
# .....\⁄ ......  .....\⁄x......      ....x.....
# .....I........  .....I........      ..I.......     
# ..............  ..............      ..........             
#                    
#                          
#  
# Em = Edge MidPoint
# I=Inner Point

LeadingPoint=[0. 0. 0.]
TrailingPoint=[0. 0. 0.]
InnerPoint   =[0. 0. 0.]    
LeadingPointOld  =[NaN NaN NaN] 
TrailingPointOld =[NaN NaN NaN] 
InnerPointOld    =[NaN NaN NaN] 
BackPoint    =[NaN NaN NaN] 

Collapsing=false
AngleSum=0.;
AreaSum=0.


#If we remove but dont clean we need to rerun the func with this flag passed out
rerunFunc=0
Idx=0
k=0
SensitiveModeOn=false
Counter=0 

#Some arrays we work on
Ang=0.
tmp=0.


NewTriP1=[0. 0. 0.]
NewTriP2=[0. 0. 0.]
NewTriP3=[0. 0. 0.]
NewTri2P1=[0. 0. 0.]
NewTri2P2=[0. 0. 0.]
NewTri2P3=[0. 0. 0.]

AngEdge_Pa=0.
AngEdge_Pb=0.
AngEdge_Pc=0.

n=Int(length(MidPoint)/3);
#Index 2 say we dont move points attached to point c around
skipindx=fill(false,n)
SplittingEdge=[0. 0. 0.;0. 0. 0.]
#Some arrays we work on
Ang=0.
tmp=0.

P1Old=copy(P1)
P2Old=copy(P2)
P3Old=copy(P3)

#For each edge loop
for i=1:length(UniqueEdges)

    b=vec(UniqueEdges[i])

    for j=b


        #Extract the points on the current bit of the edge
        (LeadingPoint,~) =GrabPointNew(LeadingPoints,P1,P2,P3,j)
        (TrailingPoint,~)=GrabPointNew(TrailingPoints,P1,P2,P3,j)
        (InnerPoint,Idx) =GrabPointNew(InnerPoints,P1,P2,P3,j)

        Pa=LeadingPoint;
        Pb=TrailingPoint;
        Pc=InnerPoint;
        PaArray=reshape(Pa,1,3)
        PbArray=reshape(Pb,1,3)
        PcArray=reshape(Pc,1,3)
        MdPnt=reshape(MidPoint[Idx,:],1,3);
        FaceNormV=reshape(FaceNormalVector[Idx,:],1,3);
        (I,FeLe,FeMd,FeEv,FeM2Ev,FeM2ELe,IntAng)=GetValues([true],PaArray,PbArray,PcArray,MdPnt,FaceNormV)


        FeEv=normr([(Pa[1]-Pb[1]) (Pa[2]-Pb[2]) (Pa[3]-Pb[3])]);
        FePa2PcV=normr([(Pc[1]-Pa[1]) (Pc[2]-Pa[2]) (Pc[3]-Pa[3])]);
        FePb2PcV=normr([(Pc[1]-Pb[1]) (Pc[2]-Pb[2]) (Pc[3]-Pb[3])]);

        #=
       #Vector po:nt:ng along edge
       FeEv=normr([(Pa[1]-Pb[1]) (Pa[2]-Pb[2]) (Pa[3]-Pb[3])]);


       FeMd=  [((Pa[1]+Pb[1])./2) ((Pa[2]+Pb[2])./2) ((Pa[3]+Pb[3])./2)];
       #Length of mid to edge dist

       FeM2ELe=sqrt.(((FeMd[1]-MidPoint[Idx,1]).^2)+((FeMd[2]-MidPoint[Idx,2]).^2)+((FeMd[3]-MidPoint[Idx,3]).^2));
       #Vector from midpoint to edge midpoint. 
       FeM2Ev=normr([FeMd[1]-MidPoint[Idx,1] FeMd[2]-MidPoint[Idx,2] FeMd[3]-MidPoint[Idx,3]]);
        =#
        

       #FeEv points from Pb to Pa 
       FeEv=FlipValue(FeEv)
           AngEdge_Pa=AngleBetweenVectors!(FePa2PcV,FeEv,tmp,AngEdge_Pa)
       FeEv=FlipValue(FeEv)

       AngEdge_Pb=AngleBetweenVectors!(FePb2PcV,FeEv,tmp,AngEdge_Pb)
       
       FePb2PcV=FlipValue(FePb2PcV);FePa2PcV=FlipValue(FePa2PcV)
           AngEdge_Pc=AngleBetweenVectors!(FePb2PcV,FePa2PcV,tmp,AngEdge_Pc)
       FePb2PcV=FlipValue(FePb2PcV);FePa2PcV=FlipValue(FePa2PcV)

       TotalAngDegrees=rad2deg(AngEdge_Pa+AngEdge_Pb+AngEdge_Pc)
       TotalAngDegrees=round(TotalAngDegrees,digits=3)
       if TotalAngDegrees!=180
           error("Internal angle is $TotalAngDegrees not 180")
       end
       
       ( Areab4,Perimb4 ) = CutAndDisplaceJulia.AreaOfTriangle3D( Pa[1],Pa[2],Pa[3],Pb[1],Pb[2],Pb[3],Pc[1],Pc[2],Pc[3] );

       NewTri1=[]
       #Now use the correct vector
       if AngEdge_Pa<AngEdge_Pb

           if rad2deg(AngEdge_Pa)<10
               println("Very low angles - bad tris")
           end
           p1=FindIntersectionOf3DVectors(FeMd,PaArray,FeM2Ev,FePa2PcV)

           #subdivide tri if new angle allows a fairly decent tri:
           NewPb2PcVec=normr([(p1[1]-Pb[1]) (p1[2]-Pb[2]) (p1[3]-Pb[3])]);
           Ang=AngleBetweenVectors!(NewPb2PcVec,FePb2PcV,tmp,Ang)
           if rad2deg(abs(Ang))>45
               NewTriP1=[NewTriP1;[p1[1] p1[2] p1[3]]]
               NewTriP2=[NewTriP2;[Pb[1] Pb[2] Pb[3]]]
               NewTriP3=[NewTriP3;[Pc[1] Pc[2] Pc[3]]]
               skipindx[Idx]=true
               #Find Connected triangle
               for k=1:3
                   NeighbourIndex=SortedTriangles[Idx,k]
                   if NeighbourIndex==0
                       continue
                   end

                   FNVTri1   = FaceNormalVector[Idx,:]  #CreateTriangleNormal(Pa,Pb,Pc)
                   FNVTri2   = FaceNormalVector[NeighbourIndex,:] #CreateTriangleNormal(Pa[NeighbourIndex,:],Pb[NeighbourIndex,:],Pc[NeighbourIndex,:])

                   #Find if it shares the edge we are interested in decimating
                   SplittingEdge[1,:].=Pc
                   SplittingEdge[2,:].=Pa 
                   InPa=ismember(SplittingEdge,P1[NeighbourIndex,:]);
                   InPb=ismember(SplittingEdge,P2[NeighbourIndex,:]);
                   InPc=ismember(SplittingEdge,P3[NeighbourIndex,:]);

                   #The edge in question
                   if sum([InPa;InPb;InPc])==2
                       #Describes if the connected edge is P1P2 P1P3 or P2P3
                       EdgeIndx=ConnectedEdge[Idx,k]
                       if EdgeIndx==12

                           #Add the new tri
                           NewTri2P1=[NewTri2P1;[p1[1] p1[2] p1[3]]]
                           NewTri2P2=[NewTri2P2;[P2[NeighbourIndex,1] P2[NeighbourIndex,2] P2[NeighbourIndex,3]]]
                           NewTri2P3=[NewTri2P3;[P3[NeighbourIndex,1] P3[NeighbourIndex,2] P3[NeighbourIndex,3]]]


                           #Change current tri in P1 P2 P3
                           #P1 the same
                           P2[NeighbourIndex,:]=[p1[1] p1[2] p1[3]]
                           #P3 the same

                       elseif EdgeIndx==13

                           #Add the new tri
                           NewTri2P1=[NewTri2P1;[p1[1] p1[2] p1[3]]]
                           NewTri2P2=[NewTri2P2;[P2[NeighbourIndex,1] P2[NeighbourIndex,2] P2[NeighbourIndex,3]]]
                           NewTri2P3=[NewTri2P3;[P3[NeighbourIndex,1] P3[NeighbourIndex,2] P3[NeighbourIndex,3]]]

                           #Change current tri in P1 P2 P3
                           #P1 the same
                           #P2 the same
                           P3[NeighbourIndex,:]=[p1[1] p1[2] p1[3]]

                       elseif EdgeIndx==23

                           #Add the new tri
                           NewTri2P1=[NewTri2P1;[P1[NeighbourIndex,1] P1[NeighbourIndex,2] P1[NeighbourIndex,3]]]
                           NewTri2P2=[NewTri2P2;[p1[1] p1[2] p1[3]]]
                           NewTri2P3=[NewTri2P3;[P3[NeighbourIndex,1] P3[NeighbourIndex,2] P3[NeighbourIndex,3]]]

                           #Change current tri in P1 P2 P3
                           #P1 the same
                           #P2 the same
                           P3[NeighbourIndex,:]=[p1[1] p1[2] p1[3]]

                       end

                       #Check normals are good
                       FNVnew=CreateTriangleNormal(NewTriP1[end,:],NewTriP2[end,:],NewTriP3[end,:])
                       Ang=AngleBetweenVectors!(FNVnew,FNVTri1,tmp,Ang)
                       if Ang>pi/2 
                           tmp=copy(NewTriP1[end,:])
                           NewTriP1[end,:]=copy(NewTriP3[end,:])
                           NewTriP3[end,:]=copy(tmp)
                       end

                       
                       FNVnew=CreateTriangleNormal(NewTri2P1[end,:],NewTri2P2[end,:],NewTri2P3[end,:])
                       Ang=AngleBetweenVectors!(FNVnew,FNVTri2,tmp,Ang)
                       if Ang>pi/2                     
                           tmp=copy(NewTri2P1[end,:])
                           NewTri2P1[end,:]=copy(NewTri2P3[end,:])
                           NewTri2P3[end,:]=copy(tmp)
                       end
                       #leave after first time 
                       #break
                   else
                       continue    
                   end

               end

           end

        (P1,P2,P3)=PushNewPoint(InnerPoints,P1,P2,P3,j,p1)
        #P3New[Idx,:]=copy(p1)

       elseif AngEdge_Pb<AngEdge_Pa

           if rad2deg(AngEdge_Pb)<10
               println("Very low angles - bad tris")
           end

           p1=FindIntersectionOf3DVectors(FeMd,PbArray,FeM2Ev,FePb2PcV)

           #subdivide tri if new angle allows a fairly decent tri:
           NewPa2PcVec=normr([(p1[1]-Pa[1]) (p1[2]-Pa[2]) (p1[3]-Pa[3])]);
           
           Ang=AngleBetweenVectors!(NewPa2PcVec,FePa2PcV,tmp,Ang)
           if rad2deg(abs(Ang))>45

               NewTriP1=[NewTriP1;[Pa[1] Pa[2] Pa[3]]]
               NewTriP2=[NewTriP2;[p1[1] p1[2] p1[3]]]
               NewTriP3=[NewTriP3;[Pc[1] Pc[2] Pc[3]]]
               skipindx[Idx]=true
               #Find Connected triangle
               for k=1:3
                   NeighbourIndex=SortedTriangles[Idx,k]
                   if NeighbourIndex==0
                       continue
                   end

                   FNVTri1   = FaceNormalVector[Idx,:] #CreateTriangleNormal(Pa,Pb,Pc)
                   FNVTri2   = FaceNormalVector[NeighbourIndex,:] #CreateTriangleNormal(Pa[NeighbourIndex,:],Pb[NeighbourIndex,:],Pc[NeighbourIndex,:])

                   #Find if it shares the edge we are interested in decimating
                   SplittingEdge[1,:].=Pc
                   SplittingEdge[2,:].=Pb
                   InPa=ismember(SplittingEdge,P1[NeighbourIndex,:]);
                   InPb=ismember(SplittingEdge,P2[NeighbourIndex,:]);
                   InPc=ismember(SplittingEdge,P3[NeighbourIndex,:]);
                   #The edge in question
                   if sum([InPa;InPb;InPc])==2
                       #Describes if the connected edge is P1P2 P1P3 or P2P3
                       EdgeIndx=ConnectedEdge[Idx,k]


                       if EdgeIndx==12

                           #Add the new tri
                           NewTri2P1=[NewTri2P1;[p1[1] p1[2] p1[3]]]
                           NewTri2P2=[NewTri2P2;[P2[NeighbourIndex,1] P2[NeighbourIndex,2] P2[NeighbourIndex,3]]]
                           NewTri2P3=[NewTri2P3;[P3[NeighbourIndex,1] P3[NeighbourIndex,2] P3[NeighbourIndex,3]]]

                           #Change current tri in P1 P2 P3
                           #P1 the same
                           P2[NeighbourIndex,:]=[p1[1] p1[2] p1[3]]
                           #P3 the same

                       elseif EdgeIndx==13

                           #Add the new tri
                           NewTri2P1=[NewTri2P1;[p1[1] p1[2] p1[3]]]
                           NewTri2P2=[NewTri2P2;[P2[NeighbourIndex,1] P2[NeighbourIndex,2] P2[NeighbourIndex,3]]]
                           NewTri2P3=[NewTri2P3;[P3[NeighbourIndex,1] P3[NeighbourIndex,2] P3[NeighbourIndex,3]]]

                           #Change current tri in P1 P2 P3
                           #P1 the same
                           #P2 the same
                           P3[NeighbourIndex,:]=[p1[1] p1[2] p1[3]]

                       elseif EdgeIndx==23

                           #Add the new tri
                           NewTri2P1=[NewTri2P1;[P1[NeighbourIndex,1] P1[NeighbourIndex,2] P1[NeighbourIndex,3]]]
                           NewTri2P2=[NewTri2P2;[p1[1] p1[2] p1[3]]]
                           NewTri2P3=[NewTri2P3;[P3[NeighbourIndex,1] P3[NeighbourIndex,2] P3[NeighbourIndex,3]]]

                           #Change current tri in P1 P2 P3
                           #P1 the same
                           #P2 the same
                           P3[NeighbourIndex,:]=[p1[1] p1[2] p1[3]]

                       end


                       #Check normals are good
                       FNVnew=CreateTriangleNormal(NewTriP1[end,:],NewTriP2[end,:],NewTriP3[end,:])

                       Ang=AngleBetweenVectors!(FNVnew,FNVTri1,tmp,Ang)
                       if Ang>pi/2 
                           tmp=copy(NewTriP1[end,:])
                           NewTriP1[end,:]=copy(NewTriP3[end,:])
                           NewTriP3[end,:]=copy(tmp)
                       end

                       
                       FNVnew=CreateTriangleNormal(NewTri2P1[end,:],NewTri2P2[end,:],NewTri2P3[end,:])
                       Ang=AngleBetweenVectors!(FNVnew,FNVTri1,tmp,Ang)
                       if Ang>pi/2                      
                           tmp=copy(NewTri2P1[end,:])
                           NewTri2P1[end,:]=copy(NewTri2P3[end,:])
                           NewTri2P3[end,:]=copy(tmp)
                       end


                       #break
                   else
                       continue    
                   end

               end

           end

           (P1,P2,P3)=PushNewPoint(InnerPoints,P1,P2,P3,j,p1)
           #P3New[Idx,:]=copy(p1)


       else


           #Use original
           #(P1,P2,P3)=PushNewPoint(InnerPoints,P1,P2,P3,j,Pc)
           #P3New[Idx,:]=copy(Pc)


       end

    end #Gone round boundary

    #Splitting certain triangles into 2 (based on area changes)
    if any(skipindx.==true)
       #Remove first part
       NewTriP1=NewTriP1[2:end,:]
       NewTriP2=NewTriP2[2:end,:]
       NewTriP3=NewTriP3[2:end,:]

       #Remove first part
       NewTri2P1=NewTri2P1[2:end,:]
       NewTri2P2=NewTri2P2[2:end,:]
       NewTri2P3=NewTri2P3[2:end,:]
    end


    for j=b


       #Extract the points on the current bit of the edge
       (InnerPoint,Idx) =GrabPointNew(InnerPoints,P1Old,P2Old,P3Old,j)    
       OldPoint=copy(InnerPoint)
       #Extract the points on the current bit of the edge
       (InnerPointNew,Idx) =GrabPointNew(InnerPoints,P1,P2,P3,j)  
       NewPoint=copy(InnerPointNew)

       if skipindx[Idx]==true
           continue
       end

       #if the connected tri contains the old point
       InP1=ismember(P1,OldPoint);
       InP2=ismember(P2,OldPoint);
       InP3=ismember(P3,OldPoint);

       #loop through and update to new point
       for k=1:length(InP1)
           if InP1[k]==true
               P1[k,:]=NewPoint;
           end
       end
       
       for k=1:length(InP2)
           if InP2[k]==true
               P2[k,:]=NewPoint;
           end
       end    
       
       for k=1:length(InP3)  
           if InP3[k]==true
               P3[k,:]=NewPoint;
           end
       end    

       #if the connected tri contains the old point
       InP1=ismember(NewTri2P1,OldPoint);
       InP2=ismember(NewTri2P2,OldPoint);
       InP3=ismember(NewTri2P3,OldPoint);

       #loop through and update to new point
       for k=1:length(InP1)
           if InP1[k]==true
               NewTri2P1[k,:]=NewPoint;
           end
       end
       
       for k=1:length(InP2)
           if InP2[k]==true
               NewTri2P2[k,:]=NewPoint;
           end
       end    
       
       for k=1:length(InP3)  
           if InP3[k]==true
               NewTri2P3[k,:]=NewPoint;
           end
       end  

   end

   #Splitting certain triangles into 2 (based on area changes)
   if any(skipindx.==true)

       P1=[P1;NewTriP1]
       P2=[P2;NewTriP2]
       P3=[P3;NewTriP3]

       P1=[P1;NewTri2P1]
       P2=[P2;NewTri2P2]
       P3=[P3;NewTri2P3]

   end

   #Resetting for next loop
   (Triangles,Points)=CreateTrianglesPointsFromP1P2P3(P1,P2,P3)
   (FaceNormalVector,MidPoint)=CreateFaceNormalAndMidPoint(Points,Triangles)
   (UniqueEdges,LeadingPoints,TrailingPoints,InnerPoints,rerunFunc,P1,P2,P3,FaceNormalVector,MidPoint)=
   GetSortedEdgesOfMeshList(P1,P2,P3,FaceNormalVector,MidPoint)
   (SortedTriangles,ConnectedEdge)=ConnectedConstraints(P1,P2,P3,MidPoint);
   n=Int(length(MidPoint)/3);
   skipindx=fill(false,n)
   P1Old=copy(P1)
   P2Old=copy(P2)
   P3Old=copy(P3)
   NewTriP1=[0. 0. 0.]
   NewTriP2=[0. 0. 0.]
   NewTriP3=[0. 0. 0.]
   NewTri2P1=[0. 0. 0.]
   NewTri2P2=[0. 0. 0.]
   NewTri2P3=[0. 0. 0.]

end

if size(P1)!=size(P2)!=size(P3)
  error("Something wrong here")
end

return P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector

end


function GrabPointNew(PointsIdxList,P1,P2,P3,j)
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


function PushNewPoint(PointsIdxList,P1,P2,P3,j,Point)
#Extract the points on the current bit of the edge

Indx=0
for k=1:3
    Idx=PointsIdxList[j,k]
    if Idx==0
        continue
    elseif k==1
        P1[Idx,:]=Point
    elseif k==2
        P2[Idx,:]=Point
    elseif k==3
        P3[Idx,:]=Point
    end
end

return P1,P2,P3
end


#=

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
    Ang[idx]=(pi/2)-(acos(dot(vec(FeIn2Ev),vec(T.FeEv))));

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
=#