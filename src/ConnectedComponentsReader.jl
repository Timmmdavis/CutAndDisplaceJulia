function ConnectedComponentsReader(Fid)

#Open the file
#using DelimitedFiles
TextFile=readdlm(Fid);

Flags=TextFile[:,1]
#Indexing of components starting from 1
Flags=Flags.+1 

Flags=Int.(Flags)

return Flags
end