function CreateTriangleNormal!(Pa,Pb,Pc,U,V,FaceNormalVector)
	
#U and V defined outside as [0 0 0]	
#FaceNormalVector as view(FaceNormalVector,i,:)

#CalculateTriangleNormal Calculates the triangle normal from 3 input
#points (edge corners of triangle). 
#Pa=(X1,Y1,Z1)
#Pb=(X2,Y2,Z2)
#Pc=(X3,Y3,Z3)

#Now calculating the normal orientation
#New vectors (see calculating normals online)
U[1]=Pb[1]-Pa[1];
U[2]=Pb[2]-Pa[2];
U[3]=Pb[3]-Pa[3];
V[1]=Pc[1]-Pa[1];
V[2]=Pc[2]-Pa[2];
V[3]=Pc[3]-Pa[3];

#Cross product of the vectors
Nx = (U[2]*V[3]) - (U[3]*V[2]);
Ny = (U[3]*V[1]) - (U[1]*V[3]);
Nz = (U[1]*V[2]) - (U[2]*V[1]);
#Vector Magnitude
aMag=sqrt((Nx * Nx) + (Ny * Ny) + (Nz * Nz));
#Norm values
FaceNormalVector[1]=Nx/aMag;
FaceNormalVector[2]=Ny/aMag;
FaceNormalVector[3]=Nz/aMag;

return FaceNormalVector
	

end