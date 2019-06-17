function ExportCrackMesh(P1,P2,P3,Dn,Dss,Dds,MeshOutputName)
#To export the actual crack shape (size opening/shearing)
#MeshOutputName - a string i.e. MeshOutputName="Mesh"
#This function guarantees that no duplicate points exist
(Triangles,Points)=CreateTrianglesPointsFromP1P2P3(P1,P2,P3)
(FaceNormalVector,MidPoint) = CreateFaceNormalAndMidPoint(Points,Triangles)
(StrikeSlipCosine,DipSlipCosine)=CalculateSSandDSDirs(FaceNormalVector[:,1],FaceNormalVector[:,2],FaceNormalVector[:,3])
(SortedTriangles,ConnectedEdge)=ConnectedConstraints(P1,P2,P3,MidPoint)

#We just add to each Point the normal vector of each connected tri, normalise this, then * by Dn
#We do the same for SS and DS bits

#=
Scl=1000;
Dn=Dn.*Scl
Dss=Dss.*Scl
Dds=Dds.*Scl
=#

#A sum of normal vectors for each point (must be normalised before using)
PointsFNVAvg=zeros(size(Points[:,2:4]))
#StrikeSlipCosine
PointsSSCAvg=zeros(size(Points[:,2:4]))
#DipSlipCosine
PointsDSCAvg=zeros(size(Points[:,2:4]))

#AverageNormal disp of each point
DnAvg=zeros(size(Points[:,1]))
#AverageNormal disp of each point
DssAvg=zeros(size(Points[:,1]))
#AverageNormal disp of each point
DdsAvg=zeros(size(Points[:,1]))

#AverageNormal disp of each point
NoConnectedTris=zeros(size(Points[:,1]))

#Compute properties over each point of the mesh
for i=1:size(Triangles,1)
	for j=1:size(Triangles,2)
		
		#Get the index of point attached to this tri
		Idx=Triangles[i,j];
		#Add the individual vectors of this tri
		PointsFNVAvg[Idx,:]+=FaceNormalVector[i,:]
		PointsSSCAvg[Idx,:]+=StrikeSlipCosine[i,:]
		PointsDSCAvg[Idx,:]+=DipSlipCosine[i,:]
		#Add the current displacements of this tri
		DnAvg[Idx]+=Dn[i]
		DssAvg[Idx]+=Dss[i]
		DdsAvg[Idx]+=Dds[i]
		#To normalise after
		NoConnectedTris[Idx]+=1

	end
end

#Drop to 0 disp for edge points of edge tris
for i=1:size(Triangles,1)
	for j=1:size(Triangles,2)
		#Get the index of point attached to this tri
		Idx=Triangles[i,j];
		#No connections to this point
		if ConnectedEdge[i,j]==0 
			#Add the current displacements of this tri
			DnAvg[Idx]=0
			DssAvg[Idx]=0
			DdsAvg[Idx]=0
		end

	end
end

#The average vectors of each point 
PointsFNVAvg=normr(PointsFNVAvg./NoConnectedTris)
PointsSSCAvg=normr(PointsSSCAvg./NoConnectedTris)
PointsDSCAvg=normr(PointsDSCAvg./NoConnectedTris)
#The average normal displacement
DnAvg=DnAvg./NoConnectedTris
DssAvg=DssAvg./NoConnectedTris
DdsAvg=DdsAvg./NoConnectedTris

#How much we displace the point by
XYZNormalDisp=PointsFNVAvg.*DnAvg
XYZStrikeSlipDisp=PointsSSCAvg.*DssAvg
XYZDipSlipDisp=PointsDSCAvg.*DdsAvg



#On the positive side of the mesh
PointsPosSide=Points[:,2:4]+XYZNormalDisp
PointsPosSide=PointsPosSide+XYZStrikeSlipDisp
PointsPosSide=PointsPosSide+XYZDipSlipDisp
#On the positive side of the mesh
PointsNegSide=Points[:,2:4]-XYZNormalDisp
PointsNegSide=PointsNegSide-XYZStrikeSlipDisp
PointsNegSide=PointsNegSide-XYZDipSlipDisp

PointsNew=[PointsPosSide;PointsNegSide]
PointsNew=[1:size(PointsNew,1) PointsNew]
TrianglesNew=[Triangles;reverse(Triangles.+maximum(Triangles),dims=2)] #reverse to flip normal on neg side

OutputDirectory=CutAndDisplaceJulia.OFFExport(PointsNew,TrianglesNew,length(TrianglesNew[:,1]),length(PointsNew[:,1]),"$MeshOutputName")

end





