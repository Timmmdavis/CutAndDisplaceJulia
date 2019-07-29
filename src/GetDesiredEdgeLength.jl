function GetDesiredEdgeLength(P1,P2,P3,NoTris)

( Area,HalfPerimeter ) = CutAndDisplaceJulia.AreaOfTriangle3D( P1[:,1],P1[:,2],P1[:,3],P2[:,1],P2[:,2],P2[:,3],P3[:,1],P3[:,2],P3[:,3] );
TotalArea=sum(Area);
DesiredTriangleArea=TotalArea/NoTris;
#Assuming equilateral tris we get the desired edge length:
target_edge_length=sqrt(DesiredTriangleArea/(sqrt(3)/4))

max_target_edge_length=maximum(HalfPerimeter)*(2/3)
#target_edge_length=(mean(HalfPerimeter)*(2/3))*1.6

return target_edge_length,max_target_edge_length
end