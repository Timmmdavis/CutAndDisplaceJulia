function TransformToADCS(y,z,by,bz,SideVec,TriVertex)

Ct=SideVec[3];
St=SideVec[2];
P1=TriVertex[2];
P2=TriVertex[3];
# Transform coordinates of the calculation points from TDCS into ADCS
(y1,z1)  =RotateObject2D(y,z,P1,P2,Ct,St)
# Transform the in-plane slip vector components from TDCS into ADCS
(by1,bz1)=RotateObject2D(by,bz,0,0,Ct,St)

return(Ct,St,y1,z1,by1,bz1)

end