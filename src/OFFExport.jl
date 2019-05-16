function OFFExport(Points,Triangles,n_tris,n_points)

#CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]))
#using Printf

#create a new file for writing
io = open("AfterAdvancingFrontClean.off", "w+");
write(io, "OFF");print(io,"\r\n");
print(io,n_points); print(io," ");print(io,n_tris); print(io," ");;print(io,0); print(io,"\r\n");
for i=1:n_points

	#13 digits precision
	@printf(io,"%0.13f", Points[i,2]);print(io," ");
	@printf(io,"%0.13f", Points[i,3]);print(io," ");
	@printf(io,"%0.13f", Points[i,4]);print(io,"\r\n");

end
for i=1:n_tris
	print(io,3);print(io," "); 
	#-1 as vertex ordering starts at 0 in .off format
	#https://people.sc.fsu.edu/~jburkardt/data/off/off.html
	print(io,Triangles[i,1]-1);print(io," ");
	print(io,Triangles[i,2]-1);print(io," ");
	print(io,Triangles[i,3]-1);print(io,"\r\n");
end
print(io,"\r\n");
close(io);

#OutputDir
OutputDirectory=string(pwd(),"\\AfterAdvancingFrontClean.off")

return OutputDirectory
end