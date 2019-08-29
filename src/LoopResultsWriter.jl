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
	@printf(io,"%e", maxX);		print(io," ");
	@printf(io,"%e", minX);		print(io," ");
	@printf(io,"%e", maxY);		print(io," ");
	@printf(io,"%e", minY);		print(io," ");
	@printf(io,"%e", maxZ);		print(io," ");
	@printf(io,"%e", minZ);		print(io," ");	

	@printf(io,"%e", G);		print(io," ");
	@printf(io,"%e", ν);		print(io," ");
	@printf(io,"%e", g);		print(io," ");
	@printf(io,"%e", Δρ);		print(io," ");
	@printf(io,"%d", NoTris);	print(io," ");
	@printf(io,"%e", KCrit);	print(io," ");	
	@printf(io,"%e", CrackVolumeIn);	print(io," ");			
	@printf(io,"%e", CrackVolumeScl);	print(io," ");			

close(io);

#OutputDir
OutputDirectory=string(pwd(),"\\$filename.txt")

return OutputDirectory
end