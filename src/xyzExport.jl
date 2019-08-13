function xyzExport(x,y,z,filename)
#vectors x y and z

#using Printf
n_points=length(x)
#create a new file for writing
io = open("$filename.xyz", "w+");
for i=1:n_points

	#13 digits precision
	@printf(io,"%0.13f", x[i]);print(io," ");
	@printf(io,"%0.13f", y[i]);print(io," ");
	@printf(io,"%0.13f", z[i]);print(io,"\r\n");

end

close(io);

if Sys.islinux()
	#OutputDir
	OutputDirectory=string(pwd(),"/$filename.xyz")
else
	#OutputDir
	OutputDirectory=string(pwd(),"\\$filename.xyz")
end

return OutputDirectory
end