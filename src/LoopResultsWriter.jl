function LoopResultsWriter(filename,PropFlag,maxX,minX,maxY,minY,maxZ,minZ,
							G,ν,g,Δρ,NoTris,KCrit,CrackVolumeIn)

#CutAndDisplaceJulia.STLExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]))
#using Printf

#create a new file for writing
io = open("$filename.txt", "w+");
write(io, "PropFlag maxX minX maxY minY maxZ minZ G ν g Δρ NoTris KCrit CrackVolumeIn");print(io,"\r\n");

	#First line
	#@printf(io,"  facet normal ")
	@printf(io,"%0.13f", PropFlag);	print(io," ");
	@printf(io,"%0.13f", maxX);		print(io," ");
	@printf(io,"%0.13f", minX);		print(io," ");
	@printf(io,"%0.13f", maxY);		print(io," ");
	@printf(io,"%0.13f", minY);		print(io," ");
	@printf(io,"%0.13f", maxZ);		print(io," ");
	@printf(io,"%0.13f", minZ);		print(io," ");	

	@printf(io,"%0.13f", G);		print(io," ");
	@printf(io,"%0.13f", ν);		print(io," ");
	@printf(io,"%0.13f", g);		print(io," ");
	@printf(io,"%0.13f", Δρ);		print(io," ");
	@printf(io,"%0.13f", NoTris);	print(io," ");
	@printf(io,"%0.13f", KCrit);	print(io," ");	
	@printf(io,"%0.13f", CrackVolumeIn);	print(io," ");			

close(io);

#OutputDir
OutputDirectory=string(pwd(),"\\$filename.txt")

return OutputDirectory
end