function IsosceliseEdgeTris(MidPoint,P1,P2,P3,Triangles,FaceNormalVector)

## PART 2: Now Rotate so always isosceles tris on edge

(FeP1P2S,FeP1P3S,FeP2P3S)=GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector);
(SortedTriangles,ConnectedEdge)=ConnectedConstraints(P1,P2,P3,MidPoint);

#Do for P1 P2: (Function at base of file)
(P1,P2,P3)=MakeEqEdgeTris(FeP1P2S,P1,P2,P3,MidPoint,FaceNormalVector,SortedTriangles,ConnectedEdge,P1,P2,P3); 

#Remesh
(Triangles,Points)=CutAndDisplaceJulia.CreateTrianglesPointsFromP1P2P3(P1,P2,P3)
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
(FeP1P2S,FeP1P3S,FeP2P3S)=GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector);
(SortedTriangles,ConnectedEdge)=ConnectedConstraints(P1,P2,P3,MidPoint);


#Do for P1 P3: (Function at base of file)
(P1,P3,P2)=MakeEqEdgeTris(FeP1P3S,P1,P3,P2,MidPoint,FaceNormalVector,SortedTriangles,ConnectedEdge,P1,P2,P3);

#Remesh
(Triangles,Points)=CutAndDisplaceJulia.CreateTrianglesPointsFromP1P2P3(P1,P2,P3)
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
(FeP1P2S,FeP1P3S,FeP2P3S)=GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector);
(SortedTriangles,ConnectedEdge)=ConnectedConstraints(P1,P2,P3,MidPoint);


#Do for P2 P3: (Function at base of file)
(P2,P3,P1)=MakeEqEdgeTris(FeP2P3S,P2,P3,P1,MidPoint,FaceNormalVector,SortedTriangles,ConnectedEdge,P1,P2,P3);

## Recreate tri
(Triangles,Points)=CreateTrianglesPointsFromP1P2P3(P1,P2,P3)
(FaceNormalVector,MidPoint) = CreateFaceNormalAndMidPoint(Points,Triangles)


return P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector

end





function MakeEqEdgeTris(T,Pa,Pb,Pc,MidPoint,FaceNormalVector,SortedTriangles,ConnectedEdge,P1,P2,P3)
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


# The aim is to isoscelise edge triangles by finding the intersection between the edge MidPointVec and 
# the longest inner edge. We then retriangulate with x as the new inner point
#
#   q                   Em            q
# ¯.¯\¯¯¯¯¯⁄¯.¯   ¯.¯\¯¯¯•¯ ⁄¯.¯      ¯.¯\¯¯ ¯¯/¯.¯
# ....\   ⁄.....  ....\  ↓⁄..... -- > ....\   /....
# .....\⁄ ....... .....\⁄x......      .....\ /.....
# .....I........  ....I........       ......x.......     
#                    
#                    
#                          
#  
# Em = Edge MidPoint
# I=Inner Point


n=Int(length(MidPoint)/3);
I=findall(T.FreeFlg);

#Vector from inner point to edge midpoint (perp 2 tri edge)
FeM2Ev=T.FeM2Ev
#Init values
FePa2PcV=Array{Float64,2}(undef, n,3);FePa2PcV=fill!(FePa2PcV, NaN)
FePb2PcV=Array{Float64,2}(undef, n,3);FePb2PcV=fill!(FePb2PcV, NaN)
FeEv=Array{Float64,2}(undef, n,3);FeEv=fill!(FeEv, NaN)

#Vector pointing along edge
FeEv[I,:]=normr([(Pa[I,1]-Pb[I,1]) (Pa[I,2]-Pb[I,2]) (Pa[I,3]-Pb[I,3])]);
FePa2PcV[I,:]=normr([(Pc[I,1]-Pa[I,1]) (Pc[I,2]-Pa[I,2]) (Pc[I,3]-Pa[I,3])]);
FePb2PcV[I,:]=normr([(Pc[I,1]-Pb[I,1]) (Pc[I,2]-Pb[I,2]) (Pc[I,3]-Pb[I,3])]);
#Get angle between edge and these vectors : smallest is the one we use
AngEdge_Pa=fill(NaN,n)
AngEdge_Pb=fill(NaN,n)
AngEdge_Pc=fill(NaN,n)
PcNew=fill(NaN,n,3)

NewTriPa=[0. 0. 0.]
NewTriPb=[0. 0. 0.]
NewTriPc=[0. 0. 0.]
NewTri2Pa=[0. 0. 0.]
NewTri2Pb=[0. 0. 0.]
NewTri2Pc=[0. 0. 0.]

#Index 2 say we dont move points attached to point c around
skipindx=fill(false,n)

for i=1:length(I)
    idx=I[i];

    #FeEv points from Pb to Pa #FeEv[I,:]=normr([(Pa[I,1]-Pb[I,1]) (Pa[I,2]-Pb[I,2]) (Pa[I,3]-Pb[I,3])]);
    AngEdge_Pa[idx]=acos(dot(vec(FePa2PcV[idx,:]),vec(.-FeEv[idx,:])))
    AngEdge_Pb[idx]=acos(dot(vec(FePb2PcV[idx,:]),vec(FeEv[idx,:])))
    AngEdge_Pc[idx]=acos(dot(vec(.-FePb2PcV[idx,:]),vec(.-FePa2PcV[idx,:])))

    TotalAngDegrees=rad2deg(AngEdge_Pa[idx]+AngEdge_Pb[idx]+AngEdge_Pc[idx])
    TotalAngDegrees=round(TotalAngDegrees,digits=3)
    if TotalAngDegrees!=180
        error("Internal angle is $TotalAngDegrees not 180")
    end
    
    ( Areab4,Perimb4 ) = CutAndDisplaceJulia.AreaOfTriangle3D( Pa[idx,1],Pa[idx,2],Pa[idx,3],Pb[idx,1],Pb[idx,2],Pb[idx,3],Pc[idx,1],Pc[idx,2],Pc[idx,3] );

    NewTri1=[]
    #Now use the correct vector
    if AngEdge_Pa[idx]<AngEdge_Pb[idx]

        if rad2deg(AngEdge_Pa[idx])<10
            println("Very low angles - bad tris")
        end
        #https://stackoverflow.com/questions/10551555/need-an-algorithm-for-3d-vectors-intersection
        
        #Func using these values
        r1=T.FeMd[idx,:] #MidPoint of the edge
        r2=Pa[idx,:]
        e1=FeM2Ev[idx,:]
        e2=FePa2PcV[idx,:];
        
        u =dot(e1,e2)
        t1=dot(r2-r1,e1)
        t2=dot(r2-r1,e2)
        d1 = (t1-u*t2)/(1.0-u*u)
       # d2 = (t2-u*t1)/(u*u.-1.0)
        p1=r1.+d1.*e1
        #p2=r2.+d2.*e2
        if u==1
            println(p1)
            println(p2)
            println("lines are parallel")
            
        end

        #subdivide tri if new angle allows a fairly decent tri:
        NewPb2PcVec=normr([(p1[1]-Pb[idx,1]) (p1[2]-Pb[idx,2]) (p1[3]-Pb[idx,3])]);
        if NewPb2PcVec!=FePb2PcV[idx,:]
            Ang=acos(dot(vec(NewPb2PcVec),vec(FePb2PcV[idx,:])))
        end
        if rad2deg(abs(Ang))>45
            NewTriPa=[NewTriPa;[Pb[idx,1] Pb[idx,2] Pb[idx,3]]]
            NewTriPb=[NewTriPb;[p1[1] p1[2] p1[3]]]
            NewTriPc=[NewTriPc;[Pc[idx,1] Pc[idx,2] Pc[idx,3]]]
            skipindx[idx]=true

            #Find Connected triangle
            for j=1:3
                NeighbourIndex=SortedTriangles[idx,j]
                if NeighbourIndex==0
                    continue
                end
                #Find if it shares the edge we are interested in decimating
                SplittingEdge=[Pc[idx,:];Pa[idx,:]] 
                InPa=ismember(SplittingEdge,Pa[NeighbourIndex,:]);
                InPb=ismember(SplittingEdge,Pb[NeighbourIndex,:]);
                InPc=ismember(SplittingEdge,Pc[NeighbourIndex,:]);
                #The edge in question
                if sum([InPa;InPb;InPc])==2
                    #Describes if the connected edge is P1P2 P1P3 or P2P3
                    EdgeIndx=ConnectedEdge[idx,j]
                    if EdgeIndx==12

                        #Add the new tri
                        NewTri2Pa=[NewTri2Pa;[P2[NeighbourIndex,1] P2[NeighbourIndex,2] P2[NeighbourIndex,3]]]
                        NewTri2Pb=[NewTri2Pb;[p1[1] p1[2] p1[3]]]
                        NewTri2Pc=[NewTri2Pc;[P3[NeighbourIndex,1] P3[NeighbourIndex,2] P3[NeighbourIndex,3]]]


                        #Change current tri in P1 P2 P3
                        #P1 the same
                        P2[NeighbourIndex,:]=[p1[1] p1[2] p1[3]]
                        #P3 the same

                    elseif EdgeIndx==13

                        #Add the new tri
                        NewTri2Pa=[NewTri2Pa;[P2[NeighbourIndex,1] P2[NeighbourIndex,2] P2[NeighbourIndex,3]]]
                        NewTri2Pb=[NewTri2Pb;[p1[1] p1[2] p1[3]]]
                        NewTri2Pc=[NewTri2Pc;[P3[NeighbourIndex,1] P3[NeighbourIndex,2] P3[NeighbourIndex,3]]]

                        #Change current tri in P1 P2 P3
                        #P1 the same
                        #P2 the same
                        P3[NeighbourIndex,:]=[p1[1] p1[2] p1[3]]

                    elseif EdgeIndx==23

                        #Add the new tri
                        NewTri2Pa=[NewTri2Pa;[P3[NeighbourIndex,1] P3[NeighbourIndex,2] P3[NeighbourIndex,3]]]
                        NewTri2Pb=[NewTri2Pb;[p1[1] p1[2] p1[3]]]
                        NewTri2Pc=[NewTri2Pc;[P1[NeighbourIndex,1] P1[NeighbourIndex,2] P1[NeighbourIndex,3]]]

                        #Change current tri in P1 P2 P3
                        #P1 the same
                        #P2 the same
                        P3[NeighbourIndex,:]=[p1[1] p1[2] p1[3]]

                    end

                    #Check normals are good
                    FNV   =CreateTriangleNormal(Pa[EdgeIndx,:],Pb[EdgeIndx,:],Pc[EdgeIndx,:])
                    FNVnew=CreateTriangleNormal(NewTriPa[end,:],NewTriPb[end,:],NewTriPc[end,:])
                    x=dot(FNVnew,FNV)
                    if abs(x)>1; Ang=0.; else; Ang=acos(x); end
                    if Ang>pi/2 
                        tmp=copy(NewTriPa[end,:])
                        NewTriPa[end,:]=copy(NewTriPc[end,:])
                        NewTriPc[end,:]=copy(tmp)
                    end

                    
                    FNV   =CreateTriangleNormal(Pa[idx,:],Pb[idx,:],Pc[idx,:])
                    FNVnew=CreateTriangleNormal(NewTri2Pa[end,:],NewTri2Pb[end,:],NewTri2Pc[end,:])
                    x=dot(FNVnew,FNV)
                    if abs(x)>1; Ang=0.; else; Ang=acos(x); end
                    if Ang>pi/2 
                        test=1
                        tmp=copy(NewTri2Pa[end,:])
                        NewTri2Pa[end,:]=copy(NewTri2Pc[end,:])
                        NewTri2Pc[end,:]=copy(tmp)
                    end

                else
                    continue    
                end

            end

        end

        PcNew[idx,:]=copy(p1)

        


    elseif AngEdge_Pb[idx]<AngEdge_Pa[idx]

        if rad2deg(AngEdge_Pb[idx])<10
            println("Very low angles - bad tris")
        end

        #Func using these values
        r1=T.FeMd[idx,:]    #MidPoint of the edge
        r2=Pb[idx,:]        #CornerPoint
        e1=FeM2Ev[idx,:]
        e2=FePb2PcV[idx,:];
        
        u =dot(e1,e2)
        t1=dot(r2-r1,e1)
        t2=dot(r2-r1,e2)
        d1 = (t1-u*t2)/(1.0-u*u)
        #d2 = (t2-u*t1)/(u*u.-1.0)
        p1=r1.+d1.*e1
        #p2=r2.+d2.*e2
        if u==1
            println(p1)
            println(p2)
            println("lines are parallel")
            
        end

        #subdivide tri if new angle allows a fairly decent tri:
        NewPa2PcVec=normr([(p1[1]-Pa[idx,1]) (p1[2]-Pa[idx,2]) (p1[3]-Pa[idx,3])]);
        if NewPa2PcVec!=FePa2PcV[idx,:]
            Ang=acos(dot(vec(NewPa2PcVec),vec(FePa2PcV[idx,:])))
        end
        if rad2deg(abs(Ang))>45
            NewTriPa=[NewTriPa;[Pa[idx,1] Pa[idx,2] Pa[idx,3]]]
            NewTriPb=[NewTriPb;[p1[1] p1[2] p1[3]]]
            NewTriPc=[NewTriPc;[Pc[idx,1] Pc[idx,2] Pc[idx,3]]]
            skipindx[idx]=true

            #Find Connected triangle
            for j=1:3
                NeighbourIndex=SortedTriangles[idx,j]
                if NeighbourIndex==0
                    continue
                end
                #Find if it shares the edge we are interested in decimating
                SplittingEdge=[Pc[idx,:];Pb[idx,:]] 
                InPa=ismember(SplittingEdge,Pa[NeighbourIndex,:]);
                InPb=ismember(SplittingEdge,Pb[NeighbourIndex,:]);
                InPc=ismember(SplittingEdge,Pc[NeighbourIndex,:]);
                #The edge in question
                if sum([InPa;InPb;InPc])==2
                    #Describes if the connected edge is P1P2 P1P3 or P2P3
                    EdgeIndx=ConnectedEdge[idx,j]
                    if EdgeIndx==12

                        #Add the new tri
                        NewTri2Pa=[NewTri2Pa;[P2[NeighbourIndex,1] P2[NeighbourIndex,2] P2[NeighbourIndex,3]]]
                        NewTri2Pb=[NewTri2Pb;[p1[1] p1[2] p1[3]]]
                        NewTri2Pc=[NewTri2Pc;[P3[NeighbourIndex,1] P3[NeighbourIndex,2] P3[NeighbourIndex,3]]]

                        #Change current tri in P1 P2 P3
                        #P1 the same
                        P2[NeighbourIndex,:]=[p1[1] p1[2] p1[3]]
                        #P3 the same

                    elseif EdgeIndx==13

                        #Add the new tri
                        NewTri2Pa=[NewTri2Pa;[P2[NeighbourIndex,1] P2[NeighbourIndex,2] P2[NeighbourIndex,3]]]
                        NewTri2Pb=[NewTri2Pb;[p1[1] p1[2] p1[3]]]
                        NewTri2Pc=[NewTri2Pc;[P3[NeighbourIndex,1] P3[NeighbourIndex,2] P3[NeighbourIndex,3]]]

                        #Change current tri in P1 P2 P3
                        #P1 the same
                        #P2 the same
                        P3[NeighbourIndex,:]=[p1[1] p1[2] p1[3]]

                    elseif EdgeIndx==23

                        #Add the new tri
                        NewTri2Pa=[NewTri2Pa;[P3[NeighbourIndex,1] P3[NeighbourIndex,2] P3[NeighbourIndex,3]]]
                        NewTri2Pb=[NewTri2Pb;[p1[1] p1[2] p1[3]]]
                        NewTri2Pc=[NewTri2Pc;[P1[NeighbourIndex,1] P1[NeighbourIndex,2] P1[NeighbourIndex,3]]]

                        #Change current tri in P1 P2 P3
                        #P1 the same
                        #P2 the same
                        P3[NeighbourIndex,:]=[p1[1] p1[2] p1[3]]

                    end

                    #Check normals are good
                    FNV   =CreateTriangleNormal(Pa[EdgeIndx,:],Pb[EdgeIndx,:],Pc[EdgeIndx,:])
                    FNVnew=CreateTriangleNormal(NewTriPa[end,:],NewTriPb[end,:],NewTriPc[end,:])
                    x=dot(FNVnew,FNV)
                    if abs(x)>1; Ang=0.; else; Ang=acos(x); end
                    if Ang>pi/2 
                        tmp=copy(NewTriPa[end,:])
                        NewTriPa[end,:]=copy(NewTriPc[end,:])
                        NewTriPc[end,:]=copy(tmp)
                    end

                    
                    FNV   =CreateTriangleNormal(Pa[idx,:],Pb[idx,:],Pc[idx,:])
                    FNVnew=CreateTriangleNormal(NewTri2Pa[end,:],NewTri2Pb[end,:],NewTri2Pc[end,:])
                    x=dot(FNVnew,FNV)
                    if abs(x)>1; Ang=0.; else; Ang=acos(x); end
                    if Ang>pi/2                      
                        tmp=copy(NewTri2Pa[end,:])
                        NewTri2Pa[end,:]=copy(NewTri2Pc[end,:])
                        NewTri2Pc[end,:]=copy(tmp)
                    end

                else
                    continue    
                end

            end

        end

        PcNew[idx,:]=copy(p1)


    else
        #Use original
        PcNew[idx,:]=copy(Pc[idx,:])
    end

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
=#

#Splitting certain triangles into 2 (based on area changes)
if any(skipindx.==true)
    #Remove first part
    NewTriPa=NewTriPa[2:end,:]
    NewTriPb=NewTriPb[2:end,:]
    NewTriPc=NewTriPc[2:end,:]

    #Remove first part
    NewTri2Pa=NewTri2Pa[2:end,:]
    NewTri2Pb=NewTri2Pb[2:end,:]
    NewTri2Pc=NewTri2Pc[2:end,:]
end

#And clean up connected tris
for i=1:length(I)

    if skipindx[I[i]]==true
        continue
    end

    OldPoint=Pc[I[i],:]
    NewPoint=PcNew[I[i],:];

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

    #if the connected tri contains the old point
    InPa=ismember(NewTri2Pa,OldPoint);
    InPb=ismember(NewTri2Pb,OldPoint);
    InPc=ismember(NewTri2Pc,OldPoint);

    #loop through and update to new point
    for j=1:length(InPa)
        if InPa[j]
            NewTri2Pa[j,:]=NewPoint;
        end
    end
    
    for j=1:length(InPb)
        if InPb[j]
            NewTri2Pb[j,:]=NewPoint;
        end
    end    
    
    for j=1:length(InPc)  
        if InPc[j]
            NewTri2Pc[j,:]=NewPoint;
        end
    end  

end

Pc[I,:]=copy(PcNew[I,:]);

#Splitting certain triangles into 2 (based on area changes)
if any(skipindx.==true)
    Pa=[Pa;NewTriPa]
    Pb=[Pb;NewTriPb]
    Pc=[Pc;NewTriPc]

    Pa=[Pa;NewTri2Pa]
    Pb=[Pb;NewTri2Pb]
    Pc=[Pc;NewTri2Pc]

end

return Pa,Pb,Pc

end

