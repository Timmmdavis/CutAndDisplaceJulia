function OFFReader(Fid)

#Open the file
#using DelimitedFiles
TextFile=readdlm(Fid);

#Number of tris
n_points=TextFile[2,1];
#Number of points
n_tris=TextFile[2,2];
#Number of edges
n_edges=TextFile[2,3];

Points=TextFile[3:2+n_points,1:3]
Triangles=TextFile[3+n_points:2+n_points+n_tris,2:4]

#Adding row with row numbers to the front of this: 
Sz=length(Points[:,1]); #getting size of rows
mono_inc=1:1:Sz; #monotomic increasing vec
Points=[mono_inc Points];

#Convert to proper type
Points=convert(Array{Float64,2},Points)
Triangles=convert(Array{Int64,2},Triangles.+1) #ordering in .off starts at 0!

return Points,Triangles
end