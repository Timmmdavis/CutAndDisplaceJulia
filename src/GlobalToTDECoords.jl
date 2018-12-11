function GlobalToTDECoords(P1,P2,P3)

(Vnorm,Vstrike,Vdip)=CalculateLocalTriCoords(P1,P2,P3)

# Transform coordinates from EFCS into TDCS
#At = [Vnorm';Vstrike';Vdip'];
Vx=[Vnorm[1],Vstrike[1],Vdip[1]];
Vy=[Vnorm[2],Vstrike[2],Vdip[2]];
Vz=[Vnorm[3],Vstrike[3],Vdip[3]];

#Init some vars
p1 = zeros(3,1);
p2 = zeros(3,1);
p3 = zeros(3,1);
P1_1=P1[1];P1_2=P1[2];P1_3=P1[3];
P2_1=P2[1];P2_2=P2[2];P2_3=P2[3];
P3_1=P3[1];P3_2=P3[2];P3_3=P3[3];

(x,y,z)=RotateObject3DNewCoords(X,Y,Z,P2_1,P2_2,P2_3,Vx,Vy,Vz)
(X1,X2,X3)=RotateObject3DNewCoords(P1_1,P1_2,P1_3,P2_1,P2_2,P2_3,Vx,Vy,Vz)
p1=[X1;X2;X3];
(X1,X2,X3)=RotateObject3DNewCoords(P3_1,P3_2,P3_3,P2_1,P2_2,P2_3,Vx,Vy,Vz)
p3=[X1;X2;X3]
p1_2=p1[2];p1_3=p1[3];
p3_2=p3[2];p3_3=p3[3];

#Get interior angles and vectors along the triangle edges. 
(e12,e13,e23,A,B,C)=CalcTDVectsAndAngles(p1,p2,p3)
	
# Determine the best arteact-free configuration for each calculation point
(casepLog,casenLog,casezLog) = trimodefinder(y,z,x,p1[2:3],p2[2:3],p3[2:3]);

return(p1,p2,p3,e12,e13,e23,A,B,C,casepLog,casenLog,casezLog,x,y,z)
end