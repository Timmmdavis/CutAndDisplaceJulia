function xyzExport(x,y,z)
#vectors x y and z

#using DelimitedFiles
writedlm("testexport.xyz",[x y z]," ")

#OutputDir
OutputDirectory=string(pwd(),"\\","testexport.xyz")

return OutputDirectory
end