function RemoveDodgyNewEdges(P1,P2,P3,NewEdgePoints,max_target_edge_length)
#After AdvancingFrontReconstruction we often end with a manifold mesh with
#large edge tris over the top and small slither tris connecting the new edge
#points. This function first removes the large triangles by assuming we only
#propagate a max of the original tri  edge length. Then it removes these small
#slithers by assuming tris that are only connected to the new edge points and
#not the original meshshould be removed

( Area,HalfPerimeter ) = AreaOfTriangle3D( P1[:,1],P1[:,2],P1[:,3],P2[:,1],P2[:,2],P2[:,3],P3[:,1],P3[:,2],P3[:,3] );
#(IntAngA,IntAngB,IntAngC)=CutAndDisplaceJulia.CalculateInternalTriAngles(P1,P2,P3)

#Clean up Advancing front result - remove overly long edges
GoodTris=HalfPerimeter.*(2/3) .< max_target_edge_length*1.5
P1=P1[GoodTris,:]
P2=P2[GoodTris,:]
P3=P3[GoodTris,:]
(Points,Triangles)=CreateTrianglesPointsFromP1P2P3(P1,P2,P3)

(FaceNormalVector,MidPoint)=CreateFaceNormalAndMidPoint(Points,Triangles)
(SortedTriangles,ConnectedEdge)=ConnectedConstraints(P1,P2,P3,MidPoint)
EdgeTri=findall(NoConnections).<3 #- edge tri

testP1=[0. 0. 0.]
testP2=[0. 0. 0.]
testP3=[0. 0. 0.]
zers  =[0. 0. 0.]
good=fill(true,size(P1,1))

#loop to remove tris that consist entirely of new edge points, we only want
#tris that are connected to the previous mesh
for i = 1:length(EdgeTri)

    ContainsNewPoint=0 #reset
    
    for j=1:length(NewEdgePoints)

        testP1.=P1[EdgeTri[i,:]].-NewEdgePoints[j,:]
        testP2.=P2[EdgeTri[i,:]].-NewEdgePoints[j,:]
        testP3.=P3[EdgeTri[i,:]].-NewEdgePoints[j,:]

        if testP1==zers || testP2==zers || testP3==zers
            ContainsNewPoint+=1
        end
        #If it only contains our new points (not the previous mesh points) its set
        #to invalid and we remove it
        if ContainsNewPoint==3
            good[EdgeTri[i]]==false
        end
    end
end

P1=copy(P1[good,:])
P2=copy(P2[good,:])
P3=copy(P3[good,:]) 

(Points,Triangles)=CreateTrianglesPointsFromP1P2P3(P1,P2,P3)

return P1,P2,P3,Points,Triangles
end