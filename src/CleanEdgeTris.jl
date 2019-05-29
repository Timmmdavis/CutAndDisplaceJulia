function CleanEdgeTris(MidPoint,P1,P2,P3,Triangles,FaceNormalVector)

#Remove slither tris
(SortedTriangles,ConnectedEdge)=ConnectedConstraints(P1,P2,P3,MidPoint)
#number of connected tris (Sorted tris rows not == to 0)
NoConnections=sum(SortedTriangles.!=0,dims=2)

#Remove any slither tris
( Area,HalfPerimeter ) = CutAndDisplaceJulia.AreaOfTriangle3D( P1[:,1],P1[:,2],P1[:,3],P2[:,1],P2[:,2],P2[:,3],P3[:,1],P3[:,2],P3[:,3] );
Good=vec(fill(true,length(Area)))
tol=mean(Area)/3
for i=1:length(Good)
    if Area[i]<tol && NoConnections[i]<2 
        Good[i]=false
    end
end

SlithersRemoved=sum(Good.==false)
println("$SlithersRemoved slither triangles have been removed")
P1=copy(P1[Good,1:3])
P2=copy(P2[Good,1:3])
P3=copy(P3[Good,1:3])
(Triangles,Points)=CreateTrianglesPointsFromP1P2P3(P1,P2,P3)
(FaceNormalVector,MidPoint) = CreateFaceNormalAndMidPoint(Points,Triangles)

rerunFunc=1 #sometime we need to run twice
extrarun=0 #for good luck
i=1
while rerunFunc==1


    (newTris,removeIndx,rerunFunc)=CutAndDisplaceJulia.CollapseEdgeTris(P1,P2,P3,MidPoint,FaceNormalVector)
    n=length(Triangles[:,1]);


    #Only doing if there are changes
    if size(removeIndx)!=()

        #Remove first dud value
        newTris=copy(newTris[2:end,:])
        removeIndx=copy(removeIndx[2:end])
        if isempty(newTris)==false
            if newTris[1,:]==newTris[end,:]
                newTris=copy(newTris[2:end,:])
            end
        end
        NoCollapsed=length(removeIndx)
        NoCreated=size(newTris,1)
        if NoCollapsed>0
            println("$NoCollapsed shared inner point edge triangles have been collapsed into $NoCreated triangles with unique inner points")
        end

        Step=collect(1:n)
        Good=vec(fill(true,n,1))
        for i=1:length(removeIndx)
            Locs=findall(in.(Step,removeIndx[i]))
            for j=1:length(Locs)
                Good[Locs[j]]=false
            end
        end

        P1=copy(P1[Good,1:3])
        P2=copy(P2[Good,1:3])
        P3=copy(P3[Good,1:3])

        P1=[P1;newTris[:,1:3]]
        P2=[P2;newTris[:,4:6]]
        P3=[P3;newTris[:,7:9]]

        ## Recreate tri
        (Triangles,Points)=CreateTrianglesPointsFromP1P2P3(P1,P2,P3)
        try (FaceNormalVector,MidPoint) = CreateFaceNormalAndMidPoint(Points,Triangles)
        catch 
            println("Check your surface, more than 2 duplicate edge tris?")
            error("Remesh here")
        end

    end
    

    #Now remove edges that have two outer edges
    (SortedTriangles,ConnectedEdge)=ConnectedConstraints(P1,P2,P3,MidPoint);
    #number of connected tris (Sorted tris rows not == to 0)
    NoConnections=sum(SortedTriangles.!=0,dims=2)
    GoodEdgeTri=NoConnections.>1 #- edge tri
    if any(GoodEdgeTri.==false)
        GoodEdgeTri=findall(vec(GoodEdgeTri)) #- edge tri
        rerunFunc=1
        P1=copy(P1[GoodEdgeTri,1:3])
        P2=copy(P2[GoodEdgeTri,1:3])
        P3=copy(P3[GoodEdgeTri,1:3])
    end

    ## Recreate tri
    (Triangles,Points)=CreateTrianglesPointsFromP1P2P3(P1,P2,P3)
    try (FaceNormalVector,MidPoint) = CreateFaceNormalAndMidPoint(Points,Triangles)
    catch 
        println("Check your surface, more than 2 duplicate edge tris?")
        error("Remesh here")
    end

    #Rerun one extra time before exit
    if rerunFunc==0
        if extrarun==0
            extrarun=1
            rerunFunc=1
        end
    end

end



return P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector
end