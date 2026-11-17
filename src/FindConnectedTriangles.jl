
# ConnectedTrianglesFinder: Finds a list of index's of connected
#                   triangles for each triangle on the surface, the
#                   distance between these and the total number of edges.
#
#                   It would be nice to optimise this as currently its a
#                   little slow. 
#               
# usage #1:
#[TotalConnectionDistance,SortedTriangles,Noconnections]...
# = ConnectedTrianglesFinder(TR,MidPoint)
#
# Arguments: (input)
#   Tr              - Output of MATLAB's triangulation function. 
#
# MidPoint          - Midpoints of the triangles of your mesh. 
#
#
# Arguments: (output)
# TotalConnectionDistance  - Each row is the distance between the triangle
#                           on this row and the other tris midpoints
#
#  SortedTriangles         - First Third and 5th column - the triangles index,
#                           other cols, the are the connected triangle
#                           index's.
#
#  Noconnections           - 1st column is the the triangles indexr, the
#                           other col the number of connected triangles
#
# Example usage (1):
#
# #If you need inputs as structured as in this BEM script use: 
#TR = triangulation(Triangles,Points(:,2:4)); 
#MidPoint= incenter(TR);
# #Call func:
#[TotalConnectionDistance,SortedTriangles,Noconnections]...
# = ConnectedTrianglesFinder(TR,MidPoint)
# # Then say you wanted to average variable StrikeSlipDisp (Dss) at every
# # tri face:
# for i=1:length(Dss)
# Mid=Dss(i); #value at tri 1
# one=Dss(SortedTriangles[i,2)); #value at 1st connected tri
# two=Dss(SortedTriangles[i,4)); #value at 2nd connected tri
# three=Dss(SortedTriangles[i,6)); #etc
# avg(i)=mean([Mid,one,two,three]);
# end
# Dss=avg;
#[ Dss] = RowVecToCol( Dss);
#
#
# Example usage (2):
# #Meshing surface to plotting the average distance
# #between tris:
#
# [x,y] = meshgrid(-2:.2:2);                                
# z = x .* exp(-x.^2 - y.^2);
# Triangles = delaunay(x(:),y(:));
# Points=[[1:length(x)]',x(:),y(:),z(:)];
# TR = triangulation(Triangles,Points(:,2:4)); 
# MidPoint= incenter(TR);
# [TotalConnectionDistance,SortedTriangles,Noconnections]...
# = ConnectedTrianglesFinder(TR,MidPoint);
# AvgConnectionDist=(sum(TotalConnectionDistance,2))./Noconnections(:,2);
# PlotSlipDistribution3d(Triangles,Points,[],AvgConnectionDist);
# 
#  Author: Tim Davis
#  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

# Open src version of : TR = triangulation(Triangles,Points(:,2),Points(:,3),Points(:,4)); E = edges(TR); %Gives all connected edges. 
# https://de.mathworks.com/matlabcentral/fileexchange/32727-fast-loop-mesh-subdivision

#(The MATLAB version was only partly translated to Julia and always stopped with an error; this
#Julia version gives the same outputs. Triangles that share an edge share two point indices.)
function FindConnectedTriangles(Triangles,MidPoint)

nNumTriangles=size(Triangles,1)

#Every edge (sorted pair of point indices) and the triangles that use it
EdgeTris=Dict{Tuple{Int,Int},Vector{Int}}()
for i=1:nNumTriangles
    for (a,b) in ((1,2),(2,3),(3,1))
        e=minmax(Int(Triangles[i,a]),Int(Triangles[i,b]))
        push!(get!(EdgeTris,e,Int[]),i)
    end
end

#Up to three connected triangles for each triangle (0 = no connection, an edge of the mesh),
#connections first
SortTri=zeros(Int,nNumTriangles,3)
for (e,tris) in EdgeTris
    length(tris)==2 || continue
    (t1,t2)=tris
    SortTri[t1,findfirst(==(0),SortTri[t1,:])]=t2
    SortTri[t2,findfirst(==(0),SortTri[t2,:])]=t1
end

#Columns 1,3,5: the triangle's index; columns 2,4,6: the connected triangles
Ascend=collect(1:nNumTriangles)
SortedTriangles=[Ascend SortTri[:,1] Ascend SortTri[:,2] Ascend SortTri[:,3]]

#Number of connections of each triangle: [index count]
Noconnections=[Ascend vec(sum(SortTri.!=0,dims=2))]

#Distance between the triangle's midpoint and each connected triangle's midpoint (0 if none)
TotalConnectionDistance=zeros(nNumTriangles,3)
for i=1:nNumTriangles, k=1:3
    j=SortTri[i,k]
    j==0 && continue
    TotalConnectionDistance[i,k]=sqrt(sum((MidPoint[i,:].-MidPoint[j,:]).^2))
end

return TotalConnectionDistance,SortedTriangles,Noconnections

end
