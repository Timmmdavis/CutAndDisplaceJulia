function STLReader(Fid)

# Fid=raw"C:\Users\user\Desktop\Sphere5120.stl"
#(Points,Triangles)=STLReader(Fid)

#Open the file
#using DelimitedFiles
TextFile=readdlm(Fid);

#Size of the array
n=length(TextFile[:,1]);

#Now create flag for vertex and triangle rows
VertexFlag=zeros(Int64, n)
VertexFlag=convert(Array{Bool,1},VertexFlag) #convert to bool
for i=1:n
	if TextFile[i,1]=="vertex"
		VertexFlag[i]=true;
	end
end	
	
	
#Extracting corresponding rows
Points=  	TextFile[VertexFlag,2:4]

#Adding row with row numbers to the front of this: 
Sz=length(Points[:,1]); #getting size of rows
mono_inc=1:1:Sz; #monotomic increasing vec
Points=[mono_inc Points];

#Creating triangles pointer
Triangles=1:length(Points[:,1]);
#Triangles=convert(Array{Float64},Triangles)
Sz2=length(Triangles)/3;
Sz2=convert(Int64,Sz2)
Triangles=reshape(Triangles,:,Sz2)';

#Convert to proper type
Points=convert(Array{Float64,2},Points)
Triangles=convert(Array{Int64,2},Triangles)

return(Points,Triangles)
end