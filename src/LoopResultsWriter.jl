function LoopResultsWriter(filename,PropFlag,maxX,minX,maxY,minY,maxZ,minZ,
							G,ν,g,Δρ,NoTris,KCrit,CrackVolumeIn,CrackVolumeScl)

#CutAndDisplaceJulia.STLExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]))
#using Printf

#@info PropFlag maxX minX maxY minY maxZ minZ G ν g Δρ NoTris KCrit CrackVolumeIn

#create a new file for writing
io = open("$filename.txt", "w+");
write(io, "PropFlag maxX minX maxY minY maxZ minZ G ν g Δρ NoTris KCrit CrackVolumeIn CrackVolumeScl");print(io,"\r\n");

	#http://www.cplusplus.com/reference/cstdio/printf/
	
	#First line
	#@printf(io,"  facet normal ")
	@printf(io,"%d", PropFlag);	print(io," ");
	@printf(io,"%f", maxX);		print(io," ");
	@printf(io,"%f", minX);		print(io," ");
	@printf(io,"%f", maxY);		print(io," ");
	@printf(io,"%f", minY);		print(io," ");
	@printf(io,"%f", maxZ);		print(io," ");
	@printf(io,"%f", minZ);		print(io," ");	

	@printf(io,"%f", G);		print(io," ");
	@printf(io,"%f", ν);		print(io," ");
	@printf(io,"%f", g);		print(io," ");
	@printf(io,"%f", Δρ);		print(io," ");
	@printf(io,"%d", NoTris);	print(io," ");
	@printf(io,"%f", KCrit);	print(io," ");	
	@printf(io,"%f", CrackVolumeIn);	print(io," ");			
	@printf(io,"%f", CrackVolumeScl);	print(io," ");			

close(io);

#OutputDir
OutputDirectory=string(pwd(),"\\$filename.txt")

return OutputDirectory
end