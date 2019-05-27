#Check we havent introduced holes by removing elements

function HoleChecker()

( Area,HalfPerimeter ) = CutAndDisplaceJulia.AreaOfTriangle3D( P1[:,1],P1[:,2],P1[:,3],P2[:,1],P2[:,2],P2[:,3],P3[:,1],P3[:,2],P3[:,3] );
(FeP1P2S,FeP1P3S,FeP2P3S)=GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector)
(SortedTriangles,ConnectedEdge)=ConnectedConstraints(P1,P2,P3,MidPoint)
#number of connected tris (Sorted tris rows not == to 0)
NoConnections=sum(SortedTriangles.!=0,dims=2)
EdgeTri=NoConnections.<3 #- edge tri

#Number of edges:
n_edges=sum([FeP1P2S.FreeFlg;FeP2P3S.FreeFlg;FeP1P3S.FreeFlg])

#FIND EDGE TRI THAT HAS NOT BEEN NULLED FROM PREVIOUS - START OF LOOP - GUARANTEES GOODNESS


#Set our flag for triangles we dont hit 
BadEdgeTri=fill(true,size(EdgeTri))
BadEdgeTri[EdgeTri==false].=false

for i=1:n_edges
	#Find next outer triangle along from this triangle (they share CurrentPoint on the outer edge)
	(triindx,CurrentPoint,TrailingPoint,CurrentEdge,InnerPoint)=
	LoopingRoundBoundary(triindx,CurrentPoint,TrailingPoint,CurrentEdge,P1,P2,P3,FeP1P2S,FeP1P3S,FeP2P3S)



end


end