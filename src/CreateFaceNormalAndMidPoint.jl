function CreateFaceNormalAndMidPoint(Points,Triangles)

#Calculate mesh surface, midpoints and normals on the triangles. 
Points2=view(Points,:,2:4);

#Prepping for loop
MidPoint=zeros(size(Triangles));
FaceNormalVector=zeros(size(Triangles));
for i = 1:length(Triangles[:,1])

	#Grabbing the current triangle points for the loop
    CurrentT1=view(Triangles,i,1);
    CurrentT2=view(Triangles,i,2);
    CurrentT3=view(Triangles,i,3);
    #Getting XYZ list for each vertex on the triangle
    Pa=view(Points2,CurrentT1,:);
    Pb=view(Points2,CurrentT2,:);
    Pc=view(Points2,CurrentT3,:);

	FaceNormalVector[i,:]=CreateTriangleNormal(Pa,Pb,Pc);
	MidPoint[i,:]=CreateMidPoint(Pa,Pb,Pc)
	
end

return FaceNormalVector,MidPoint

end