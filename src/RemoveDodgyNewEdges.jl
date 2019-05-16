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
P1=copy(P1[GoodTris,:])
P2=copy(P2[GoodTris,:])
P3=copy(P3[GoodTris,:])
NoRemoved=sum(GoodTris.==false)
println("Removed $NoRemoved faces from AdvancingFront as their area was too large")



testP1=[0. 0. 0.]
testP2=[0. 0. 0.]
testP3=[0. 0. 0.]
zers  =[0. 0. 0.]
good=fill(true,size(P1,1))


#loop to remove tris that consist entirely of new edge points, we only want
#tris that are connected to the previous mesh
for i = 1:size(P1,1)

    ContainsNewPoint=0 #reset
    
    for j=1:size(NewEdgePoints,1)

        testP1=P1[i,:].-NewEdgePoints[j,:]
        testP2=P2[i,:].-NewEdgePoints[j,:]
        testP3=P3[i,:].-NewEdgePoints[j,:]

        if testP1==zers
            @bp
            ContainsNewPoint+=1
        end 

        if testP2==zers 
            @bp
            ContainsNewPoint+=1
        end
        if testP3==zers
            @bp
            ContainsNewPoint+=1
            println("Inside")
        end


            
        #If it only contains our new points (not the previous mesh points) its set
        #to invalid and we remove it
        if ContainsNewPoint==3
            println("SuperGood")
            good[EdgeTri[i]]==false
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
