function xyzExport(x,y,z,filename)
#vectors x y and z

#using DelimitedFiles
writedlm("$filename.xyz",[x y z]," ")

#OutputDir
OutputDirectory=string(pwd(),"\\",filename,".xyz")

return OutputDirectory
end