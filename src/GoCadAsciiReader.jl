function GoCadAsciiReader(Fid)

# Fid=raw"C:\Users\user\Desktop\SphereUniformDistributionON_46Faces.ts"
#(Points,Triangles)=GoCadAsciiReader(Fid)

#Open the file
#using DelimitedFiles
TextFile=readdlm(Fid);

#Size of the array
n=length(TextFile[:,1]);

#Check first col doesnt contain the string "PROJECTION"
for i=1:n
	if TextFile[i,1]=="PROJECTION"
		 error("More than one Gocad ascii surface in the file, this function can only deal with 1")
	end
end

#Now create flag for vertex and triangle rows
VertexFlag=zeros(Int64, n)
TriangleFlag=zeros(Int64, n)
AtomFlag=zeros(Int64, n)
VertexFlag=convert(Array{Bool,1},VertexFlag) #convert to bool
TriangleFlag=convert(Array{Bool,1},TriangleFlag) #convert to bool
AtomFlag=convert(Array{Bool,1},AtomFlag) #convert to bool
for i=1:n
	if TextFile[i,1]=="VRTX"
		VertexFlag[i]=true;
	end
	if TextFile[i,1]=="TRGL"
		TriangleFlag[i]=true;
	end	
	if TextFile[i,1]=="ATOM"
		AtomFlag[i]=true;
	end		
end	
	
	
#Extracting corresponding rows
Points=  	TextFile[VertexFlag,2:5]
Triangles=  TextFile[TriangleFlag,2:4]
Atoms=      TextFile[AtomFlag,2:3]

#Now if atoms exist we need to plug these in to the points vector:
if any(Atoms);
    Tag=zeros(length(Atoms));
    Atoms=[Atoms Tag]; #Adding rows to be filled
    #Looping and filling in correct XYZ for the atom numbers
    for i=1:length(Atoms[:,1])
        #Find where points row 1 matches atom num
        a=Points[:,1]==Atoms[i,2];
        #Push in correct row
        Atoms[i,2:4]=Points[a,2:4];   
    end
    #Now appending atoms onto points and sorting this by row 1
    Points=[Points;Atoms];
	Points=sortslices(Points; dims=1);
end    


#Fixing a bug, this codebase relies on the fact the cols are the same as
#numbers the triangles are pointing to in the first row of 'Points'.
#Some program exports do not get this right. 
for i=1:length(Points[:,1])
    if Points[i,1] != i          #if the row is not equal to the first column on that row. 
        #any numbers correponding to this - 'Points(i,1)' need to be changed
        #to this 'i' in triangles
        NumToChange=Points[i,1];
        Bad = Triangles == NumToChange;   
        Triangles[Bad]=i;
        Points[i,1]=i;
    end
end   
#Convert to proper type
Points=convert(Array{Float64,2},Points)
Triangles=convert(Array{Int64,2},Triangles)

return(Points,Triangles)

end