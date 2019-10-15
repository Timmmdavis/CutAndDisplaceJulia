function PutOnEdgeOfCirc(MidPoint,P1,P2,P3,FaceNormalVector,Radius,X,Y)
#Assuming its a penny shaped mesh we make sure all the edge points lie on the edge of a real circle.
#X and Y are the index's of the XY coords, so even if its vertical you just supply these as 1 and 3


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
lps=0
#For each edge loop
for i=1:length(UniqueEdges)

  #Run through list of this edge (Unique edges is a unit range)
  #Add 10 more loops to unit range to conver first tri - This assumes we dont
  #have a crazy edge tri with loads of small segments. Less than 10 if that loop only has < 10 segments
  maxlength=length(UniqueEdges[i])-1
  if maxlength>10
    maxlength=10
  end
  b=[vec(UniqueEdges[i]);vec(minimum(UniqueEdges[i]):minimum(UniqueEdges[i])+maxlength)]

  for j=b
    #Get the two triangles surrounding j
    j_last=j-1

    if j_last==0
      continue
    end

    j_now=j
    j_future=j+1
    if j==1
      j_trailing=b[end]
    elseif j==b[end]
      j_future=b[1]
    end

    #Extract the points on the current bit of the edge
    (k_now,Idx_trailing) =GrabPointNew7(TrailingPoints,P1,P2,P3,j_now)
    (k_last,Idx_leading) =GrabPointNew7(LeadingPoints,P1,P2,P3,j_last)

    if k_now==1
      CurrentPoint=P1[Idx_trailing,:]
    elseif k_now==2
      CurrentPoint=P2[Idx_trailing,:]
    elseif k_now==3
      CurrentPoint=P3[Idx_trailing,:]
    end     

    if k_last==1
      CurrentPoint2=P1[Idx_leading,:]
    elseif k_last==2
      CurrentPoint2=P2[Idx_leading,:]
    elseif k_last==3
      CurrentPoint2=P3[Idx_leading,:]
    end   


    if CurrentPoint==CurrentPoint2


      if k_now==1
        P1Old=copy(P1[Idx_trailing,:])
        P1=ApplyConstraints(P1,Idx_trailing,X,Y,Radius)
        #println(P1[Idx_trailing,:])

        for i=1:3
          if SortedTriangles[Idx_trailing,i]>0
            IndxOfConn=SortedTriangles[Idx_trailing,i]
            if ConnectedEdge[Idx_trailing,i]==12
              if P1[IndxOfConn,:]==P1Old
                P1[IndxOfConn,:]=P1[Idx_trailing,:]
              elseif P2[IndxOfConn,:]==P1Old
                P2[IndxOfConn,:]=P1[Idx_trailing,:]
              end
            elseif ConnectedEdge[Idx_trailing,i]==13
              if P1[IndxOfConn,:]==P1Old
                P1[IndxOfConn,:]=P1[Idx_trailing,:]
              elseif P3[IndxOfConn,:]==P1Old
                P3[IndxOfConn,:]=P1[Idx_trailing,:]
              end
            elseif ConnectedEdge[Idx_trailing,i]==23
              if P2[IndxOfConn,:]==P1Old
                P2[IndxOfConn,:]=P1[Idx_trailing,:]
              elseif P3[IndxOfConn,:]==P1Old
                P3[IndxOfConn,:]=P1[Idx_trailing,:]
              end
            end
          end
        end

      elseif k_now==2
        P2Old=copy(P2[Idx_trailing,:])
        P2=ApplyConstraints(P2,Idx_trailing,X,Y,Radius)
        #println(P1[Idx_trailing,:])

        for i=1:3
          if SortedTriangles[Idx_trailing,i]>0
            IndxOfConn=SortedTriangles[Idx_trailing,i]
            if ConnectedEdge[Idx_trailing,i]==12
              if P1[IndxOfConn,:]==P2Old
                P1[IndxOfConn,:]=P2[Idx_trailing,:]
              elseif P2[IndxOfConn,:]==P2Old
                P2[IndxOfConn,:]=P2[Idx_trailing,:]
              end
            elseif ConnectedEdge[Idx_trailing,i]==13
              if P1[IndxOfConn,:]==P2Old
                P1[IndxOfConn,:]=P2[Idx_trailing,:]
              elseif P3[IndxOfConn,:]==P2Old
                P3[IndxOfConn,:]=P2[Idx_trailing,:]
              end
            elseif ConnectedEdge[Idx_trailing,i]==23
              if P2[IndxOfConn,:]==P2Old
                P2[IndxOfConn,:]=P2[Idx_trailing,:]
              elseif P3[IndxOfConn,:]==P2Old
                P3[IndxOfConn,:]=P2[Idx_trailing,:]
              end
            end
          end
        end

      elseif k_now==3
        P3Old=copy(P3[Idx_trailing,:])
        P3=ApplyConstraints(P3,Idx_trailing,X,Y,Radius)
        #println(P1[Idx_trailing,:])

        for i=1:3
          if SortedTriangles[Idx_trailing,i]>0
            IndxOfConn=SortedTriangles[Idx_trailing,i]
            if ConnectedEdge[Idx_trailing,i]==12
              if P1[IndxOfConn,:]==P3Old
                P1[IndxOfConn,:]=P3[Idx_trailing,:]
              elseif P2[IndxOfConn,:]==P3Old
                P2[IndxOfConn,:]=P3[Idx_trailing,:]
              end
            elseif ConnectedEdge[Idx_trailing,i]==13
              if P1[IndxOfConn,:]==P3Old
                P1[IndxOfConn,:]=P3[Idx_trailing,:]
              elseif P3[IndxOfConn,:]==P3Old
                P3[IndxOfConn,:]=P3[Idx_trailing,:]
              end
            elseif ConnectedEdge[Idx_trailing,i]==23
              if P2[IndxOfConn,:]==P3Old
                P2[IndxOfConn,:]=P3[Idx_trailing,:]
              elseif P3[IndxOfConn,:]==P3Old
                P3[IndxOfConn,:]=P3[Idx_trailing,:]
              end
            end
          end
        end
      end     

      
      if k_last==1
        P1Old=copy(P1[Idx_leading,:])
        P1=ApplyConstraints(P1,Idx_leading,X,Y,Radius)
        #println(P1[Idx_leading,:])
        for i=1:3
          if SortedTriangles[Idx_leading,i]>0
            IndxOfConn=SortedTriangles[Idx_leading,i]
            if ConnectedEdge[Idx_leading,i]==12
              if P1[IndxOfConn,:]==P1Old
                P1[IndxOfConn,:]=P1[Idx_leading,:]
              elseif P2[IndxOfConn,:]==P1Old
                P2[IndxOfConn,:]=P1[Idx_leading,:]
              end
            elseif ConnectedEdge[Idx_leading,i]==13
              if P1[IndxOfConn,:]==P1Old
                P1[IndxOfConn,:]=P1[Idx_leading,:]
              elseif P3[IndxOfConn,:]==P1Old
                P3[IndxOfConn,:]=P1[Idx_leading,:]
              end
            elseif ConnectedEdge[Idx_leading,i]==23
              if P2[IndxOfConn,:]==P1Old
                P2[IndxOfConn,:]=P1[Idx_leading,:]
              elseif P3[IndxOfConn,:]==P1Old
                P3[IndxOfConn,:]=P1[Idx_leading,:]
              end
            end
          end
        end
      elseif k_last==2
        P2Old=copy(P2[Idx_leading,:])
        P2=ApplyConstraints(P2,Idx_leading,X,Y,Radius)
        #println(P2[Idx_leading,:])
        for i=1:3
          if SortedTriangles[Idx_leading,i]>0
            IndxOfConn=SortedTriangles[Idx_leading,i]
            if ConnectedEdge[Idx_leading,i]==12
              if P1[IndxOfConn,:]==P2Old
                P1[IndxOfConn,:]=P2[Idx_leading,:]
              elseif P2[IndxOfConn,:]==P2Old
                P2[IndxOfConn,:]=P2[Idx_leading,:]
              end
            elseif ConnectedEdge[Idx_leading,i]==13
              if P1[IndxOfConn,:]==P2Old
                P1[IndxOfConn,:]=P2[Idx_leading,:]
              elseif P3[IndxOfConn,:]==P2Old
                P3[IndxOfConn,:]=P2[Idx_leading,:]
              end
            elseif ConnectedEdge[Idx_leading,i]==23
              if P2[IndxOfConn,:]==P2Old
                P2[IndxOfConn,:]=P2[Idx_leading,:]
              elseif P3[IndxOfConn,:]==P2Old
                P3[IndxOfConn,:]=P2[Idx_leading,:]
              end
            end
          end
        end
      elseif k_last==3
        P3Old=copy(P3[Idx_leading,:])
        P3=ApplyConstraints(P3,Idx_leading,X,Y,Radius)
        #println(P3[Idx_leading,:])
        for i=1:3
          if SortedTriangles[Idx_leading,i]>0
            IndxOfConn=SortedTriangles[Idx_leading,i]
            if ConnectedEdge[Idx_leading,i]==12
              if P1[IndxOfConn,:]==P3Old
                P1[IndxOfConn,:]=P3[Idx_leading,:]
              elseif P2[IndxOfConn,:]==P3Old
                P2[IndxOfConn,:]=P3[Idx_leading,:]
              end
            elseif ConnectedEdge[Idx_leading,i]==13
              if P1[IndxOfConn,:]==P3Old
                P1[IndxOfConn,:]=P3[Idx_leading,:]
              elseif P3[IndxOfConn,:]==P3Old
                P3[IndxOfConn,:]=P3[Idx_leading,:]
              end
            elseif ConnectedEdge[Idx_leading,i]==23
              if P2[IndxOfConn,:]==P3Old
                P2[IndxOfConn,:]=P3[Idx_leading,:]
              elseif P3[IndxOfConn,:]==P3Old
                P3[IndxOfConn,:]=P3[Idx_leading,:]
              end
            end
          end
        end  
      end   
      

    end


  end

end 

(Triangles,Points)=CutAndDisplaceJulia.CreateTrianglesPointsFromP1P2P3(P1,P2,P3)
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )

return P1,P2,P3,MidPoint,FaceNormalVector,Points,Triangles

end

function GrabPointNew7(PointsIdxList,P1,P2,P3,j)
#Extract the points on the current bit of the edge
InnerPointNo=0;
Indx=0
for k=1:3
    Idx=PointsIdxList[j,k]
    if Idx==0
        continue
    elseif k==1
        InnerPointNo=k
        Indx=Idx
    elseif k==2
        InnerPointNo=k
        Indx=Idx
    elseif k==3
        InnerPointNo=k
        Indx=Idx
    end
end

return InnerPointNo,Indx
end


function ApplyConstraints(Px,i,X,Y,Radius)

OldPoint=copy(Px[i,:])

x=Px[i,X];
y=Px[i,Y];
(θ,ρ)=CutAndDisplaceJulia.cart2pol(x,y)
(x,y)=CutAndDisplaceJulia.pol2cart(θ,Radius)
Px[i,X]=x;
Px[i,Y]=y;

return Px

end