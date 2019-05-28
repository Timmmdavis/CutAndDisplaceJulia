function CollapseEdgeTris2(P1,P2,P3,MidPoint,FaceNormalVector)
## The aim is to collapse edge triangles that share an inner point
# We do this by looping round and checking if tris share a point inside the surface
#													Ev
#   T   C 	 	       T   C 		          T   ---> C
# ¯.¯\¯¯|¯¯/¯.¯  ¯.¯\¯¯|¯¯/¯.¯       ¯.¯\¯¯ ¯¯/¯.¯
# ....\ | /....  ....\ | /....  -- > ....\   /....
# .....\|/.....  .....\|/.....       .....\ /.....
# ......I......  ......I......       ......I......   
#	    
#  
# T=TrailingPoint
# C=CurrentPoint (or LeadingPoint)
# I=InnerPoint
# Ev=EdgeVector

(UniqueEdges,LeadingPoints,TrailingPoints,InnerPoints)=
GetSortedEdgesOfMeshList(P1,P2,P3,FaceNormalVector,MidPoint)

( Area,HalfPerimeter ) = CutAndDisplaceJulia.AreaOfTriangle3D( P1[:,1],P1[:,2],P1[:,3],P2[:,1],P2[:,2],P2[:,3],P3[:,1],P3[:,2],P3[:,3] );

LeadingPoint=[0. 0. 0.]
TrailingPoint=[0. 0. 0.]
InnerPoint	 =[0. 0. 0.]	
LeadingPointOld  =[NaN NaN NaN]	
TrailingPointOld =[NaN NaN NaN]	
InnerPointOld    =[NaN NaN NaN]	
BackPoint    =[NaN NaN NaN]	

Collapsing=false
AngleSum=0.;
AreaSum=0.

#Index's of triangles we will remove after loop
removeIndx=-1000; #To be removed
#Memoryless vector to fill with the list of new points
fillme=zeros(1,9)
#Appending to this vector everytime we add a new tri
newTris=fill(NaN,1,9)
#If we remove but dont clean we need to rerun the func with this flag passed out
rerunFunc=0
Idx=0

@bp

#For each edge loop
for i=1:length(UniqueEdges)
	#Run through list of this edge (Unique edges is a unit range)
	for j=UniqueEdges[i] 

		LeadingPointOld  =LeadingPoint
		TrailingPointOld =TrailingPoint
		InnerPointOld    =InnerPoint
		IdxOld=Idx

		#Extract the points on the current bit of the edge
		(LeadingPoint,~) =GrabPoint(LeadingPoints,P1,P2,P3,j)
		(TrailingPoint,~)=GrabPoint(TrailingPoints,P1,P2,P3,j)
		(InnerPoint,Idx) =GrabPoint(InnerPoints,P1,P2,P3,j)
		
		@bp 
		
		#Check if this edge shares the inner point with the old one -if so we collapse it
		if InnerPoint==InnerPointOld

			#Set the back point and reset some vars if we have just started collapsing
			if Collapsing==false
				BackPoint=TrailingPointOld
				#Set the indx to remove inside
				if IdxOld!=0
					removeIndx=[removeIndx IdxOld] #add to list of tris to remove
				end
				#reset some values
				V1=normr([LeadingPoint[1]-InnerPoint[1]  LeadingPoint[2]-InnerPoint[2]  LeadingPoint[3]-InnerPoint[3]])
				V2=normr([TrailingPoint[1]-InnerPoint[1] TrailingPoint[2]-InnerPoint[2] TrailingPoint[3]-InnerPoint[3]])
				x=dot(V1,V2)
				if abs(x)>1; Ang=0.; else; Ang=acos(x); end					  			
				AngleSum=Ang
				AreaSum=0.;
			end
			#
			removeIndx=[removeIndx Idx] #add to list of tris to remove

			###############################################################
			#Vectors pointing out from inner point on current tri
			V1=normr([LeadingPoint[1]-InnerPoint[1]  LeadingPoint[2]-InnerPoint[2]  LeadingPoint[3]-InnerPoint[3]])
			V2=normr([TrailingPoint[1]-InnerPoint[1] TrailingPoint[2]-InnerPoint[2] TrailingPoint[3]-InnerPoint[3]])
			x=dot(V1,V2)
			if abs(x)>1; Ang=0.; else; Ang=acos(x); end					  			
			AngleSum+=Ang
			#Sum of areas of each connected tri
			AreaSum+=Area[Idx]
			###############################################################

			#We are now collapsing the tri
			Collapsing=true

		#The latest inner point doesnt match so create the new tri and start again
		elseif InnerPoint!=InnerPointOld && Collapsing==true

			Collapsing=false #reset to false

			#If we get here we know that there is no edge connection to the current triangle
			#Adding to list of new triangles
			fillme[1:3].=InnerPoint;
			fillme[4:6].=LeadingPoint;
			fillme[7:9].=BackPoint;

			#make sure point ordering matches that of the normal direction
			AvgVect = FaceNormalVector[Idx,:]
			NewTriNormalVector = CreateTriangleNormal( fillme[1:3],fillme[4:6],fillme[7:9] );
			AllignFlag = dot(AvgVect,NewTriNormalVector);
			if AllignFlag<0
				#flip order
			    fillme[1:3].=BackPoint;
			    fillme[7:9].=InnerPoint;
			end  

			##Check that the area is retained, if not the new triangle is just bad
			( NewArea,~ ) = CutAndDisplaceJulia.AreaOfTriangle3D( fillme[:,1],fillme[:,2],fillme[:,3],fillme[:,4],fillme[:,5],fillme[:,6],fillme[:,7],fillme[:,8],fillme[:,9] );
			#AngleSum - if we have passed through over 180 degrees its bound to be bad
			if NewArea[1]<(AreaSum/4) || rad2deg(AngleSum)>180 #way too small
				@info NewArea[1] (AreaSum/4) rad2deg(AngleSum)
				#Do nothing 
				rerunFunc=1 #we will need to rerun
				#break
			else
				println("Whoop")
				newTris=[newTris; copy(fillme) ]
			end

		end


			
	end
end

return newTris,removeIndx,rerunFunc

end


function GrabPoint(PointsIdxList,P1,P2,P3,j)
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
