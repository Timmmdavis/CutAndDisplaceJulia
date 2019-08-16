#Get the lastest png in the file:

#cd(raw"C:\Users\timmm\Desktop\MeshProp\WhereTheMeshesLive")


function FindingLastPng(Dir)
x=readdir(Dir)

beststring=0;
MaxCounter=0;
for i=1:length(x)
	FileExtension=split(x[i],'.')
	if FileExtension[end]=="png"
		Counter=split(FileExtension[1],'-')
		Counter=parse(Int,Counter[1])
		if Counter>MaxCounter
			beststring=i;
			MaxCounter=Counter
		end
	end
end

return(x[beststring])
end