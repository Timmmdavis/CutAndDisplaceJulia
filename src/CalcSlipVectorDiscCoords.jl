function CalcSlipVectorDiscCoords(SideVec,eZ,X,Y,Z,PA,bX,bY,bZ,beta)
#Calculate the Slip vector in dislocation coordinates

ey1 = [SideVec[1:2];0];
ey1 = ey1/norm(ey1);
ey3 = -eZ;
ey2 = cross(ey3,ey1);

# Transform coordinates from EFCS to the first ADCS
(y1A,y2A,y3A)=RotateObject3DNewCoords(X,Y,Z,PA[1],PA[2],PA[3],ey1,ey2,ey3)

# Transform coordinates from EFCS to the second ADCS
(y1AB,y2AB,y3AB)=RotateObject3DNewCoords(SideVec[1],SideVec[2],SideVec[3],0,0,0,ey1,ey2,ey3)
y1B = y1A.-y1AB;
y2B = y2A.-y2AB;
y3B = y3A.-y3AB;

# Transform slip vector components from EFCS to ADCS
(b1,b2,b3)=RotateObject3DNewCoords(bX,bY,bZ,0,0,0,ey1,ey2,ey3)

# Determine the best arteact-free configuration for the calculation
# points near the free furface
I = (beta.*y1A).>=0;

return(b1,b2,b3,I,y1A,y2A,y3A,y1B,y2B,y3B,ey1,ey2,ey3)
end
