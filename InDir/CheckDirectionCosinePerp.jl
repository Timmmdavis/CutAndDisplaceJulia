function CheckDirectionCosinePerp(C1,C2,C3)
#Checks the direction cosines are all perpendicular, If not, something is wrong. 
#C1=DipSlipCosine;
#C2=StrikeSlipCosine;
#C3=FaceNormalVector;

#Extracting parts of the new coordinates. 
A11=C1[:,1];A12=C1[:,2];A13=C1[:,3];
A21=C2[:,1];A22=C2[:,2];A23=C2[:,3];
A31=C3[:,1];A32=C3[:,2];A33=C3[:,3];

DotPro12=A11.*A12+A21.*A22+A31.*A32;
DotPro13=A11.*A13+A21.*A23+A31.*A33;
DotPro23=A12.*A13+A22.*A23+A32.*A33;

@info DotPro12 DotPro13 DotPro23
if DotPro12[1]>eps() || DotPro13[1]>eps() || DotPro23[1]>eps() 
	error("Direction cosines are not at 90 degrees to each other")
end    

end