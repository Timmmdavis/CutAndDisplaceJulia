function LoadData(ModuleName,DataNameAsString)
# LoadData Creates a string that is the path to a file within a modules test directory.
# Works with Windows and Linux

#Load triangles and points from file (mesh)
ModuleDir=pathof(ModuleName);
ModuleDir=splitdir(ModuleDir); #remove file name
ModuleDir=ModuleDir[1];
ModuleDir=splitdir(ModuleDir); #out of src
ModuleDir=ModuleDir[1];
if Sys.iswindows()
    SurfaceDir=string(ModuleDir,string("\\test\\",DataNameAsString))
else
	SurfaceDir=string(ModuleDir,string("/test/",DataNameAsString))
end


end