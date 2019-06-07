function CreateFaceNormalAndMidPoint(Points,Triangles)

#Calculate mesh surface, midpoints and normals on the triangles. 
Points2=view(Points,:,2:4);

U=[0.; 0.; 0.]
V=[0.; 0.; 0.]
Pa=[0.; 0.; 0.]
Pb=[0.; 0.; 0.]
Pc=[0.; 0.; 0.]
FNV=[0.; 0.; 0.]
MP=[0.; 0.; 0.]
CurrentT1=[0]
CurrentT2=[0]
CurrentT3=[0]


#Prepping for loop
MidPoint=zeros(size(Triangles));
FaceNormalVector=zeros(size(Triangles));

n=length(Triangles[:,1]);

for i = 1:n

	#Grabbing the current triangle points for the loop
    CurrentT1[1]=Triangles[i,1];
    CurrentT2[1]=Triangles[i,2];
    CurrentT3[1]=Triangles[i,3];
	
    #Getting XYZ list for each vertex on the triangle
    Pa[1]=Points2[CurrentT1[1],1];
    Pa[2]=Points2[CurrentT1[1],2];
    Pa[3]=Points2[CurrentT1[1],3];
    Pb[1]=Points2[CurrentT2[1],1];
    Pb[2]=Points2[CurrentT2[1],2];
    Pb[3]=Points2[CurrentT2[1],3];
    Pc[1]=Points2[CurrentT3[1],1];
    Pc[2]=Points2[CurrentT3[1],2];
    Pc[3]=Points2[CurrentT3[1],3];

    FNV[1]=FaceNormalVector[i,1];
    FNV[2]=FaceNormalVector[i,2];
    FNV[3]=FaceNormalVector[i,3];
    MP[1]=MidPoint[i,1];
    MP[2]=MidPoint[i,2];
    MP[3]=MidPoint[i,3];

    FNV=CreateTriangleNormal!(Pa,Pb,Pc,U,V,FNV);
    MP=CreateMidPoint!(Pa,Pb,Pc,MP)
    
    FaceNormalVector[i,1]=FNV[1]
    FaceNormalVector[i,2]=FNV[2]
    FaceNormalVector[i,3]=FNV[3]
    MidPoint[i,1]=MP[1]
    MidPoint[i,2]=MP[2]
    MidPoint[i,3]=MP[3]
    

end

return FaceNormalVector,MidPoint

end

function CreateFaceNormalAndMidPoint(P1::Array,P2::Array,P3::Array)

U=[0.; 0.; 0.]
V=[0.; 0.; 0.]
Pa=[0.; 0.; 0.]
Pb=[0.; 0.; 0.]
Pc=[0.; 0.; 0.]
FNV=[0.; 0.; 0.]
MP=[0.; 0.; 0.]

#Prepping for loop
MidPoint=zeros(size(P1));
FaceNormalVector=zeros(size(P1));
n=length(P1[:,1]);

for i = 1:n

    #Getting XYZ list for each vertex on the triangle
    Pa[1]=P1[i,1];
    Pa[2]=P1[i,2];
    Pa[3]=P1[i,3];
    Pb[1]=P2[i,1];
    Pb[2]=P2[i,2];
    Pb[3]=P2[i,3];
    Pc[1]=P3[i,1];
    Pc[2]=P3[i,2];
    Pc[3]=P3[i,3];

    FNV[1]=FaceNormalVector[i,1];
    FNV[2]=FaceNormalVector[i,2];
    FNV[3]=FaceNormalVector[i,3];
    MP[1]=MidPoint[i,1];
    MP[2]=MidPoint[i,2];
    MP[3]=MidPoint[i,3];

    FNV=CreateTriangleNormal!(Pa,Pb,Pc,U,V,FNV);
    MP=CreateMidPoint!(Pa,Pb,Pc,MP)

    FaceNormalVector[i,1]=FNV[1]
    FaceNormalVector[i,2]=FNV[2]
    FaceNormalVector[i,3]=FNV[3]
    MidPoint[i,1]=MP[1]
    MidPoint[i,2]=MP[2]
    MidPoint[i,3]=MP[3]
    
end

return FaceNormalVector,MidPoint

end