function RemoveDodgyNewEdges(P1,P2,P3,Points,Triangles,FaceNormalVector,MidPoint,OnlyEdgePoints,max_target_edge_length)
#After AdvancingFrontReconstruction we often end with a manifold mesh with
#large edge tris over the top and small slither tris connecting the new edge
#points. This function first removes the large triangles by assuming we only
#propagate a max of the original tri  edge length. Then it removes these small
#slithers by assuming tris that are only connected to the new edge points and
#not the original meshshould be removed

#=
( Area,HalfPerimeter ) = AreaOfTriangle3D( P1[:,1],P1[:,2],P1[:,3],P2[:,1],P2[:,2],P2[:,3],P3[:,1],P3[:,2],P3[:,3] );
#(IntAngA,IntAngB,IntAngC)=CutAndDisplaceJulia.CalculateInternalTriAngles(P1,P2,P3)

#Clean up Advancing front result - remove overly long edges
GoodTris=HalfPerimeter.*(2/3) .< max_target_edge_length*1.5
P1=copy(P1[GoodTris,:])
P2=copy(P2[GoodTris,:])
P3=copy(P3[GoodTris,:])
NoRemoved=sum(GoodTris.==false)
println("Removed $NoRemoved faces from AdvancingFront as their area was too large")
=#

(SortedTriangles,ConnectedEdge)=ConnectedConstraints(P1,P2,P3,MidPoint)
#number of connected tris (Sorted tris rows not == to 0)
NoConnections=sum(SortedTriangles.!=0,dims=2)
EdgeTri=NoConnections.<3 #- edge tri

good=fill(true,size(P1,1))

#Only use unique ones
OnlyEdgePoints=unique(OnlyEdgePoints,dims=1)
#All my exports to .off and .xyz etc are in 13 trailing point digits
OnlyEdgePoints=round.(OnlyEdgePoints,digits=13)
P1=round.(P1,digits=13)
P2=round.(P2,digits=13)
P3=round.(P3,digits=13)



#loop to remove tris that consist entirely of only edge points, we only want
#tris that are connected to the previous mesh
for i = 1:size(P1,1)

    NoEdgePoints=0 #reset
    
    for j=1:size(OnlyEdgePoints,1)

        testP1=P1[i,:]==OnlyEdgePoints[j,:]
        testP2=P2[i,:]==OnlyEdgePoints[j,:]
        testP3=P3[i,:]==OnlyEdgePoints[j,:]

        if testP1
            NoEdgePoints+=1
        end 
        if testP2
            NoEdgePoints+=1
        end
        if testP3
            NoEdgePoints+=1
        end
    end
    #If it only contains our new points (not the previous mesh points) its set
    #to invalid and we remove it
    if NoEdgePoints==3
        good[i]=false
        continue
    end

    #Part 2 - remove edges with over 90 degrees between connected face normal
    #vectors We only do when more bad connectons than good. The propagation
    #algorithm doesnt allow for over 90 degree turns. This is because in rare
    #cases the algorithm to remove bad tris above doesnt quite catch all the
    #baddies. Sometimes new tris go from outer edge new points to old outer
    #tris in a dodgy way
    if NoEdgePoints>0
        if EdgeTri[i] #if its an edge triangle
            
            GoodCons=0; #reset
            BadCons=0;

            for j=1:NoConnections[i]

                ConnectedNormal=FaceNormalVector[SortedTris[i,j],:]
                CurrentNormal=FaceNormalVector[i,:]
                AngBetweenNormals=acos(dot(vec(ConnectedNormal),vec(CurrentNormal)))
                
                if AngBetweenNormals>(pi/2) #bigger than 90
                    BadCons+=1;
                else    #smaller
                    GoodCons+=1
                end

            end    
            #Only remove if its mainly got bad connections
            if BadCons>GoodCons 
                good[i]=false
            end
        end    
    end

end



P1=copy(P1[good,:])
P2=copy(P2[good,:])
P3=copy(P3[good,:]) 
NoRemoved=sum(good.==false)
println("Removed $NoRemoved faces from AdvancingFront as they were only the new tip points")

(Points,Triangles)=CreateTrianglesPointsFromP1P2P3(P1,P2,P3)

return P1,P2,P3,Points,Triangles
end
