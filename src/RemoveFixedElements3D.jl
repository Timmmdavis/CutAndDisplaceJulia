function RemoveFixedElements3D(FixedEls,Triangles,FaceNormalVector,MidPoint,P1,P2,P3)

#Get logical flag of bits to keep
FD=FixedEls.==0;
#Needs to be a single vector
FD=reshape(FD,length(FD))
#Remove bits. 
Triangles=Triangles[FD,:];
FaceNormalVector=FaceNormalVector[FD,:];
MidPoint=MidPoint[FD,:];
P1=P1[FD,:];
P2=P2[FD,:];
P3=P3[FD,:];

return FixedEls,Triangles,FaceNormalVector,MidPoint,P1,P2,P3

end