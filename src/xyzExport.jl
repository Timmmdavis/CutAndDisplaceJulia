function xyzExport(Points)

#using DelimitedFiles
writedlm("testexport.txt",Points[:,2:4]," ")

#OutputDir
OutputDirectory=pwd()

return OutputDirectory
end