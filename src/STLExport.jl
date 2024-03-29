function STLExport(P1,P2,P3,FaceNormalVector,filename)

#CutAndDisplaceJulia.STLExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]))
#using Printf

#create a new file for writing
io = open("$filename.stl", "w+");
write(io, "solid STL generated by Herr Davis");print(io,"\r\n");
for i=1:size(P1,1)

	#First line
	@printf(io,"  facet normal ")
	@printf(io,"%0.13f", FaceNormalVector[i,1]);print(io," ");
	@printf(io,"%0.13f", FaceNormalVector[i,2]);print(io," ");
	@printf(io,"%0.13f", FaceNormalVector[i,3]);
	print(io,"\r\n");

	#NextLine
	@printf(io,"    outer loop");
	print(io,"\r\n");
	
	#NextLine
	@printf(io,"      vertex ");
	@printf(io,"%0.13f", P1[i,1]);print(io," ");
	@printf(io,"%0.13f", P1[i,2]);print(io," ");
	@printf(io,"%0.13f", P1[i,3]);
	print(io,"\r\n");	
	#NextLine
	@printf(io,"      vertex ");
	@printf(io,"%0.13f", P2[i,1]);print(io," ");
	@printf(io,"%0.13f", P2[i,2]);print(io," ");
	@printf(io,"%0.13f", P2[i,3]);
	print(io,"\r\n");	
	#NextLine
	@printf(io,"      vertex ");
	@printf(io,"%0.13f", P3[i,1]);print(io," ");
	@printf(io,"%0.13f", P3[i,2]);print(io," ");
	@printf(io,"%0.13f", P3[i,3]);
	print(io,"\r\n");	

	#NextLine
	@printf(io,"    endloop");
	print(io,"\r\n");

	#NextLine
	@printf(io,"  endfacet");
	print(io,"\r\n");

end

@printf(io,"endsolid vcg");
print(io,"\r\n");
close(io);

if Sys.islinux()
	#OutputDir
	OutputDirectory=string(pwd(),"/$filename.stl")
else
	#OutputDir
	OutputDirectory=string(pwd(),"\\$filename.stl")
end



return OutputDirectory
end