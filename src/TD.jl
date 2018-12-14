function TD(X,Y,Z,P1List,P2List,P3List,Dss,Dds,Dn,nu,mu,DispFlag,StrainFlag,HSflag)
# TDdispFS 
# Calculates displacements associated with a triangular dislocation in an 
# elastic full-space.
#
# TD: Triangular Dislocation
# EFCS: Earth-Fixed Coordinate System
# TDCS: Triangular Dislocation Coordinate System
# ADCS: Angular Dislocation Coordinate System
# 
# INPUTS
# X, Y and Z: 
# Coordinates of calculation points in EFCS (East, North, Up). X, Y and Z 
# must have the same size.
#
# P1,P2 and P3:
# Coordinates of TD vertices in EFCS.
# 
# Ss, Ds and Ts:
# TD slip vector components (Strike-slip, Dip-slip, Tensile-slip).
#
# nu:
# Poisson's ratio.
#
# OUTPUTS
# Ux, Uy and Uz:
# Calculated displacement vector components in EFCS. Ux, Uy and Uz have
# the same Uyit as Ss, Ds and Ts in the inputs.
# 
# 
# Example: Calculate and plot the first component of displacement vector 
# on a regular grid.
# 
# (X,Y,Z) = meshgrid(-3:.02:3,-3:.02:3,2);
# (Ux,Uy,Uz) = TDdispFS(X,Y,Z,(-1 0 0),(1 -1 -1),(0 1.5 .5),-1,2,3,.25);
# h = surf(X,Y,reshape(Ux,size(X)),'edgecolor','none');
# view(2)
# axis equal
# axis tight
# set(gcf,'renderer','painters')
#
# Reference journal article: 
# Nikkhoo M. and Walter T.R., 2015. Triangular dislocation: An analytical, 
# artefact-free solution. 
# Submitted to Geophysical Journal International 
#
# Copyright (c) 2014 Mehdi Nikkhoo
# 
# Permission is hereby granted, free of charge, to any person obtaining a 
# copy of this software and associated documentation files 
# (the "Software"), to deal in the Software without restriction, including 
# without limitation the rights to use, copy, modify, merge, publish, 
# distribute, sublicense, and/or sell copies of the Software, and to permit
# persons to whom the Software is furnished to do so, subject to the 
# following conditions:
# 
# The above copyright notice and this permission notice shall be included 
# in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
# NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
# USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# I appreciate any comments or bug reports.
#
# Mehdi Nikkhoo
# created: 2012.5.14
# Last modified: 2014.7.30
# 
# VolcanoTectonics Research Group
# Section 2.1, Physics of Earthquakes and Volcanoes
# Department 2, Physics of the Earth
# Helmholtz Centre Potsdam
# German Research Centre for Geosciences (GFZ)
# 
# email: 
# mehdi.nikkhoo@gfz-potsdam.de 
# mehdi.nikkhoo@gmail.com



SzCmp=size(P1List,1); if length(P1List)==3; SzCmp=1; end

CorrectDimsFlg= SzCmp==size(P2List,1) &&
				SzCmp==size(P3List,1) &&
				SzCmp==length(Dds) &&
				SzCmp==length(Dss) &&
				SzCmp==length(Dn);
if CorrectDimsFlg!=1
	error("Element size inputs must be the same dimensions")
end



println("This initial allocation should follow through the functions....")
#Do some allocation before loop
if DispFlag==1
	UxDn = Array{Float64}(undef, length(X),SzCmp); 
	UyDn = Array{Float64}(undef, length(X),SzCmp); 
	UzDn = Array{Float64}(undef, length(X),SzCmp); 
	UxDss = Array{Float64}(undef, length(X),SzCmp); 
	UyDss = Array{Float64}(undef, length(X),SzCmp); 
	UzDss = Array{Float64}(undef, length(X),SzCmp); 
	UxDds = Array{Float64}(undef, length(X),SzCmp); 
	UyDds = Array{Float64}(undef, length(X),SzCmp); 
	UzDds = Array{Float64}(undef, length(X),SzCmp); 
end
if StrainFlag==1
	ExxDn = Array{Float64}(undef, length(X),SzCmp); 
	EyyDn = Array{Float64}(undef, length(X),SzCmp); 
	EzzDn = Array{Float64}(undef, length(X),SzCmp);
	ExyDn = Array{Float64}(undef, length(X),SzCmp); 
	ExzDn = Array{Float64}(undef, length(X),SzCmp); 
	EyzDn = Array{Float64}(undef, length(X),SzCmp);
	ExxDss = Array{Float64}(undef, length(X),SzCmp); 
	EyyDss = Array{Float64}(undef, length(X),SzCmp); 
	EzzDss = Array{Float64}(undef, length(X),SzCmp);
	ExyDss = Array{Float64}(undef, length(X),SzCmp); 
	ExzDss = Array{Float64}(undef, length(X),SzCmp); 
	EyzDss = Array{Float64}(undef, length(X),SzCmp);
	ExxDds = Array{Float64}(undef, length(X),SzCmp); 
	EyyDds = Array{Float64}(undef, length(X),SzCmp); 
	EzzDds = Array{Float64}(undef, length(X),SzCmp);
	ExyDds = Array{Float64}(undef, length(X),SzCmp); 
	ExzDds = Array{Float64}(undef, length(X),SzCmp); 
	EyzDds = Array{Float64}(undef, length(X),SzCmp);
end

for i=1:SzCmp #For every element

	println("Remove this line once strain done")
	bx = Dn[i];  # DisplacementNormal
	by = Dss[i]; # DisplacementStrike-slip
	bz = Dds[i]; # DisplacementDip-slip
	P1=P1List[i,:];
	P2=P2List[i,:];
	P3=P3List[i,:];

	if HSflag==1
		if any(Z .>0) | any(P1[3] .>0) | any(P2[3] .>0) | any(P3[3] .>0)
			error("Half-space solution: Z coordinates must be negative!")
		end
		P1i=copy(P1);
		P2i=copy(P2);
		P3i=copy(P3);
		P1i[3] = -P1[3];
		P2i[3] = -P2[3];
		P3i[3] = -P3[3];
	end

	if DispFlag==1

		# Calculate main dislocation contribution to displacements
		(UxDn[:,i],UyDn[:,i],UzDn[:,i],
		 UxDss[:,i],UyDss[:,i],UzDss[:,i],
		 UxDds[:,i],UyDds[:,i],UzDds[:,i]) = 
		 TDdispFS(X,Y,Z,P1,P2,P3,Dss[i],Dds[i],Dn[i],nu);
		
		if HSflag==1
			
			
			# Calculate harmonic fUyction contribution to displacements
			(UxDnFSC,UyDnFSC,UzDnFSC,
			 UxDssFSC,UyDssFSC,UzDssFSC,
			 UxDdsFSC,UyDdsFSC,UzDdsFSC) = 
			 TDdisp_HarFunc(X,Y,Z,P1,P2,P3,Dss[i],Dds[i],Dn[i],nu);

			# Calculate image dislocation contribution to displacements
			(UxDnIS,UyDnIS,UzDnIS,
			 UxDssIS,UyDssIS,UzDssIS,
			 UxDdsIS,UyDdsIS,UzDdsIS) =
			 TDdispFS(X,Y,Z,P1i,P2i,P3i,Dss[i],Dds[i],Dn[i],nu);
			if P1i[3]==0 && P2i[3]==0 && P3i[3]==0
				UzDnIS  = -UzDnIS;
				UzDssIS = -UzDssIS;
				UzDdsIS = -UzDdsIS;
			end

			
			# Calculate the complete displacement vector components in EFCS
			UxDn[:,i]  = UxDn[:,i]+UxDnIS+UxDnFSC;
			UyDn[:,i]  = UyDn[:,i]+UyDnIS+UyDnFSC;
			UzDn[:,i]  = UzDn[:,i]+UzDnIS+UzDnFSC;
			UxDss[:,i] = UxDss[:,i]+UxDssIS+UxDssFSC;
			UyDss[:,i] = UyDss[:,i]+UyDssIS+UyDssFSC;
			UzDss[:,i] = UzDss[:,i]+UzDssIS+UzDssFSC;
			UxDds[:,i] = UxDds[:,i]+UxDdsIS+UxDdsFSC;
			UyDds[:,i] = UyDds[:,i]+UyDdsIS+UyDdsFSC;
			UzDds[:,i] = UzDds[:,i]+UzDdsIS+UzDdsFSC;
			
			if P1i[3]==0 && P2i[3]==0 && P3i[3]==0
				UxDn[:,i]  = -UxDn[:,i];
				UyDn[:,i]  = -UyDn[:,i];
				UzDn[:,i]  = -UzDn[:,i];
				UxDss[:,i] = -UxDss[:,i];
				UyDss[:,i] = -UyDss[:,i];
				UzDss[:,i] = -UzDss[:,i];
				UxDds[:,i] = -UxDds[:,i];
				UyDds[:,i] = -UyDds[:,i];
				UzDds[:,i] = -UzDds[:,i];				
			end
			
		end
		
	end 

	if StrainFlag==1

		#Elastic con
		lambda=(2*mu*nu)/(1-(2*nu));
		
		(ExxDn[:,i],EyyDn[:,i],EzzDn[:,i],ExyDn[:,i],ExzDn[:,i],EyzDn[:,i],
		 ExxDss[:,i],EyyDss[:,i],EzzDss[:,i],ExyDss[:,i],ExzDss[:,i],EyzDss[:,i],
		 ExxDds[:,i],EyyDds[:,i],EzzDds[:,i],ExyDds[:,i],ExzDds[:,i],EyzDds[:,i]) =
		 TDstrainFS(X,Y,Z,P1,P2,P3,Dss[i],Dds[i],Dn[i],mu,lambda)	
		
		if HSflag==1
		
			(ExxFSC,EyyFSC,EzzFSC,ExyFSC,ExzFSC,EyzFSC) = TDstrain_HarFunc(X,Y,Z,P1,P2,P3,by,bz,bx,mu,lambda,nu);
		
			# Calculate image dislocation contribution to strains and stresses
			(ExxIS,EyyIS,EzzIS,ExyIS,ExzIS,EyzIS) = TDstrainFS(X,Y,Z,P1i,P2i,P3i,by,bz,bx,mu,lambda);
			if P1i[3]==0 && P2i[3]==0 && P3i[3]==0
				ExzIS = -ExzIS;
				EyzIS = -EyzIS;
			end
			
			Exx[:,i]=Exx[:,i]+ExxFSC+ExxIS;
			Eyy[:,i]=Eyy[:,i]+EyyFSC+EyyIS;
			Ezz[:,i]=Ezz[:,i]+EzzFSC+EzzIS;
			Exy[:,i]=Exy[:,i]+ExyFSC+ExyIS;
			Exz[:,i]=Exz[:,i]+ExzFSC+ExzIS;
			Eyz[:,i]=Eyz[:,i]+EyzFSC+EyzIS;
			
		end #HS flag
		
	end #Stress flag
	
end #Over all els
	
	
if DispFlag==0
	return(Exx,Eyy,Ezz,Exy,Exz,Eyz)
elseif StrainFlag==0
	return(UxDn,UyDn,UzDn,
		   UxDss,UyDss,UzDss,
		   UxDds,UyDds,UzDds)
else

	return(ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn,
		   ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
		   ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,
		   UxDn,UyDn,UzDn,
		   UxDss,UyDss,UzDss,
		   UxDds,UyDds,UzDds)
	
end
end


function TDdispFS(X,Y,Z,P1,P2,P3,Dss,Dds,Dn,nu)
#Calculate displacement components in a full space
(Vnorm,Vstrike,Vdip)=CalculateLocalTriCoords(P1,P2,P3)
(p1,p2,p3,e12,e13,e23,A,B,C,casepLog,casenLog,casezLog,x,y,z)=GlobalToTDECoords(P1,P2,P3,X,Y,Z,Vnorm,Vstrike,Vdip)
p1_2=p1[2];p1_3=p1[3];
p3_2=p3[2];p3_3=p3[3];

xn=x[casenLog];
yn=y[casenLog];
zn=z[casenLog];
xp=x[casepLog];
yp=y[casepLog];
zp=z[casepLog];

# Calculate first angular dislocation contribution POS
(ux1Tp,uy1Tp,uz1Tp,vx1Tp,vy1Tp,vz1Tp,wx1Tp,wy1Tp,wz1Tp) = TDSetupD(xp,yp,zp,A,Dn,Dss,Dds,nu,p1,-e13);
# Calculate second angular dislocation contribution
(ux2Tp,uy2Tp,uz2Tp,vx2Tp,vy2Tp,vz2Tp,wx2Tp,wy2Tp,wz2Tp) = TDSetupD(xp,yp,zp,B,Dn,Dss,Dds,nu,p2,e12);
# Calculate third angular dislocation contribution
(ux3Tp,uy3Tp,uz3Tp,vx3Tp,vy3Tp,vz3Tp,wx3Tp,wy3Tp,wz3Tp) = TDSetupD(xp,yp,zp,C,Dn,Dss,Dds,nu,p3,e23);	

# Calculate first angular dislocation contribution NEG
(ux1Tn,uy1Tn,uz1Tn,vx1Tn,vy1Tn,vz1Tn,wx1Tn,wy1Tn,wz1Tn) = TDSetupD(xn,yn,zn,A,Dn,Dss,Dds,nu,p1,e13);
# Calculate second angular dislocation contribution
(ux2Tn,uy2Tn,uz2Tn,vx2Tn,vy2Tn,vz2Tn,wx2Tn,wy2Tn,wz2Tn) = TDSetupD(xn,yn,zn,B,Dn,Dss,Dds,nu,p2,-e12);
# Calculate third angular dislocation contribution
(ux3Tn,uy3Tn,uz3Tn,vx3Tn,vy3Tn,vz3Tn,wx3Tn,wy3Tn,wz3Tn) = TDSetupD(xn,yn,zn,C,Dn,Dss,Dds,nu,p3,-e23);		

#Do some allocation before loop
println("Send these directly into the func above? Add in there (reduce alloc massiv amounts)")
uxV = Array{Float64}(undef, length(x),1); 
vxV = Array{Float64}(undef, length(x),1); 
wxV = Array{Float64}(undef, length(x),1);
uyV = Array{Float64}(undef, length(x),1); 
vyV = Array{Float64}(undef, length(x),1); 
wyV = Array{Float64}(undef, length(x),1);
uzV = Array{Float64}(undef, length(x),1); 
vzV = Array{Float64}(undef, length(x),1); 
wzV = Array{Float64}(undef, length(x),1);

uxV[casenLog]=ux1Tn+ux2Tn+ux3Tn;
vxV[casenLog]=vx1Tn+vx2Tn+vx3Tn;
wxV[casenLog]=wx1Tn+wx2Tn+wx3Tn;
uxV[casepLog]=ux1Tp+ux2Tp+ux3Tp;
vxV[casepLog]=vx1Tp+vx2Tp+vx3Tp;
wxV[casepLog]=wx1Tp+wx2Tp+wx3Tp;

uyV[casenLog]=uy1Tn+uy2Tn+uy3Tn;
vyV[casenLog]=vy1Tn+vy2Tn+vy3Tn;
wyV[casenLog]=wy1Tn+wy2Tn+wy3Tn;
uyV[casepLog]=uy1Tp+uy2Tp+uy3Tp;
vyV[casepLog]=vy1Tp+vy2Tp+vy3Tp;
wyV[casepLog]=wy1Tp+wy2Tp+wy3Tp;

uzV[casenLog]=uz1Tn+uz2Tn+uz3Tn;
vzV[casenLog]=vz1Tn+vz2Tn+vz3Tn;
wzV[casenLog]=wz1Tn+wz2Tn+wz3Tn;
uzV[casepLog]=uz1Tp+uz2Tp+uz3Tp;
vzV[casepLog]=vz1Tp+vz2Tp+vz3Tp;
wzV[casepLog]=wz1Tp+wz2Tp+wz3Tp;

# Calculate the "incomplete" displacement vector components in TDCS
for i=1:length(x)

	if casezLog[i] == 1; 
		uxV[i] = NaN;
		vxV[i] = NaN;
		wxV[i] = NaN;
		uyV[i] = NaN;
		vyV[i] = NaN;
		wyV[i] = NaN;
		uzV[i] = NaN;
		vzV[i] = NaN;
		wzV[i] = NaN;
	end

	a1=-x[i];	a2=p1_2[1]-y[i];	a3=p1_3[1]-z[i];
	b1=-x[i];	b2=-y[i];		b3=-z[i];
	c1=-x[i];	c2=p3_2[1]-y[i];	c3=p3_3[1]-z[i];
	a12=a1^2;
	# Calculate the Burgers' function contribution corresponding to the TD
	na = sqrt((a12+a2.^2+a3.^2));
	nb = sqrt((a12+b2.^2+b3.^2));
	nc = sqrt((a12+c2.^2+c3.^2));
	one=a1*(b2*c3-b3*c2)-
		a2*(b1*c3-b3*c1)+
		a3*(b1*c2-b2*c1);
	ab=a1*b1+a2*b2+a3*b3;
	ac=a1*c1+a2*c2+a3*c3;
	bc=b1*c1+b2*c2+b3*c3;
	two=na*nb*nc+ab*nc+ac*nb+bc*na;
	Fi = -2*atan(one,two)/4/pi; 
	
	# Calculate the complete displacement vector components in TDCS
	uxV[i] = Dn*Fi+uxV[i];
	vyV[i] = Dss*Fi+vyV[i];
	wzV[i] = Dds*Fi+wzV[i];
	#Only add parts that matter

	
end	

println("Remove allocation of new vars below")
# Transform the complete displacement vector components from TDCS into EFCS
(UxDn,UyDn,UzDn)=   RotateObject3DNewCoords(uxV,vxV,wxV,0,0,0,Vnorm,Vstrike,Vdip)
(UxDss,UyDss,UzDss)=RotateObject3DNewCoords(uyV,vyV,wyV,0,0,0,Vnorm,Vstrike,Vdip)
(UxDds,UyDds,UzDds)=RotateObject3DNewCoords(uzV,vzV,wzV,0,0,0,Vnorm,Vstrike,Vdip)

return(UxDn,UyDn,UzDn,UxDss,UyDss,UzDss,UxDds,UyDds,UzDds)
end

function GlobalToTDECoords(P1,P2,P3,X,Y,Z,Vnorm,Vstrike,Vdip)
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
(p1[1],p1[2],p1[3])=RotateObject3DNewCoords(P1_1,P1_2,P1_3,P2_1,P2_2,P2_3,Vx,Vy,Vz)
(p3[1],p3[2],p3[3])=RotateObject3DNewCoords(P3_1,P3_2,P3_3,P2_1,P2_2,P2_3,Vx,Vy,Vz)#Vx,Vy,Vz
p1_2=p1[2];p1_3=p1[3];
p3_2=p3[2];p3_3=p3[3];

#Get interior angles and vectors along the triangle edges. 
(e12,e13,e23,A,B,C)=CalcTDVectsAndAngles(p1,p2,p3)
	
# Determine the best arteact-free configuration for each calculation point
(casepLog,casenLog,casezLog) = trimodefinder(y,z,x,p1[2:3],p2[2:3],p3[2:3]);

return(p1,p2,p3,e12,e13,e23,A,B,C,casepLog,casenLog,casezLog,x,y,z)
end

function CalculateLocalTriCoords(P1,P2,P3)
# Calculate unit strike, dip and normal to TD vectors: For a horizontal TD 
# as an exception, if the normal vector points upward, the strike and dip 
# vectors point Northward and Westward, whereas if the normal vector points
# downward, the strike and dip vectors point Southward and Westward, 
# respectively.

eY = [0;1;0];
eZ = [0;0;1];

Vnorm = cross(P2-P1,P3-P1);
Vnorm = Vnorm/norm(Vnorm);
Vstrike = cross(eZ,Vnorm);
if norm(Vstrike)==0
	Vstrike = eY*Vnorm[3];
end
Vstrike = Vstrike/norm(Vstrike);
Vdip = cross(Vnorm,Vstrike);

return(Vnorm,Vstrike,Vdip)

end

function CalcTDVectsAndAngles(p1,p2,p3)
#Get interior angles and vectors along the triangle edges. 
#Inputs:
#p1,p2,p3 	- 1*3 vectors containing xyz of each triangles corner point. 
#Outputs
#e12 etc 	- unit vectors along the TD sides in the triangles coordinate system (norm strike dip)
#Outputs
#A B C 		- TD angles interior angles. 

# Calculate the unit vectors along TD sides in TDCS
e12=p2-p1;
e12=e12/norm(e12);
e13=p3-p1;
e13=e13/norm(e13);
e23=p3-p2;
e23=e23/norm(e23);

# Calculate the TD angles
A=e12'*e13;
A=acos(A[1]);
B=-e12'*e23;
B=acos(B[1]);
C=e23'*e13;
C=acos(C[1]);

return(e12,e13,e23,A,B,C)
end

function trimodefinder(x,y,z,p1,p2,p3)
# trimodefinder calculates the normalized barycentric coordinates of 
# the points with respect to the TD vertices and specifies the appropriate
# artefact-free configuration of the angular dislocations for the 
# calculations. The input matrices x, y and z share the same size and
# correspond to the y, z and x coordinates in the TDCS, respectively. p1,
# p2 and p3 are two-component matrices representing the y and z coordinates
# of the TD vertices in the TDCS, respectively.
# The components of the output (trimode) corresponding to each calculation 
# points, are 1 for the first configuration, -1 for the second 
# configuration and 0 for the calculation point that lie on the TD sides.

#Get index's
p1_1=p1[1];
p2_1=p2[1];
p3_1=p3[1];
p1_2=p1[2];
p2_2=p2[2];
p3_2=p3[2];
#Init some values outside of the loop
p2_2_m_p3_2=p2_2-p3_2;
p3_1_m_p2_1=p3_1-p2_1;
p1_1_m_p3_1=p1_1-p3_1;
p1_2_m_p3_2=p1_2-p3_2;
p3_2_m_p1_2=p3_2-p1_2;
Base=(p2_2_m_p3_2*p1_1_m_p3_1+p3_1_m_p2_1*p1_2_m_p3_2);

trimode=ones(Int64, length(x),1)
for i=1:length(x)

	a = (p2_2_m_p3_2*(x[i]-p3_1)+p3_1_m_p2_1*(y[i]-p3_2))/Base;
	b = (p3_2_m_p1_2*(x[i]-p3_1)+p1_1_m_p3_1*(y[i]-p3_2))/Base;
	c = 1-a-b;
		
	if a<=0 && b>c && c>a
		trimode[i] = -1;
	elseif b<=0 && c>a && a>b
		trimode[i] = -1;
	elseif c<=0 && a>b && b>c	
		trimode[i] = -1;
	elseif a==0 && b>=0 && c>=0		
		trimode[i] = 0;
	elseif a>=0 && b==0 && c>=0
		trimode[i] = 0;
	elseif a>=0 && b>=0 && c==0
		trimode[i] = 0;	
	end
		
	if trimode[i]==0 && z!=0
		trimode[i] = 1;
	end

end
casepLog=falses(length(trimode),1);
casenLog=falses(length(trimode),1);
casezLog=falses(length(trimode),1);
for i=1:length(trimode)
	if trimode[i]==1
		casepLog[i] = true; 
	end
	if trimode[i]==-1
		casenLog[i] = true; 
	end	
	if trimode[i]==0;
		casezLog[i] = true;
	end
end

return(casepLog,casenLog,casezLog)
end


function TDSetupD(x,y,z,alpha,Dn,Dss,Dds,nu,TriVertex,SideVec)
# TDSetupD transforms coordinates of the calculation points as well as 
# slip vector components from ADCS into TDCS. It then calculates the 
# displacements in ADCS and transforms them into TDCS.

(Ct,St,y1,z1,Dss1,Dds0)=TransformToADCS(y,z,1.,0.,SideVec,TriVertex)
(Dss0,Dds1)=RotateObject2D(0.,1.,0,0,Ct,St)

#Init arrays
Ang=-pi+alpha;
cosA = cos(Ang);
sinA = sin(Ang);

println("Wont need to be allocated if passed directly in")
Ux = Array{Float64}(undef, length(x),1);
Uy = Array{Float64}(undef, length(x),1);
Uz = Array{Float64}(undef, length(x),1);
Vx = Array{Float64}(undef, length(x),1);
Vy = Array{Float64}(undef, length(x),1);
Vz = Array{Float64}(undef, length(x),1);
Wx = Array{Float64}(undef, length(x),1);
Wy = Array{Float64}(undef, length(x),1);
Wz = Array{Float64}(undef, length(x),1);

#Extra defs out of loop to speed it up
E1=(1-nu); #Elastic cons
E2=(1-2*nu);
cosA2=cosA^2;
sinADE1=sinA/8/pi/(1-nu);

Dn8p=Dn/8/pi;

# Calculate displacements associated with an angular dislocation in ADCS
for i=1:length(x)
	
	(ux,uy,uz,vx,vy,vz,wx,wy,wz) = AngDisDisp(x[i],y1[i],z1[i],cosA,sinA,E1,E2,cosA2,sinADE1);
	
	#components due to opening
	Ux[i]=Dn8p/E1*ux;
	Vx[i]=Dn8p/E1*vx;
	Wx[i]=Dn8p/E1*wx;
	#comp mixed components (in the current coords)
	Uy[i]=(Dss1/8/pi/E1*uy)+(Dds0*sinADE1*uz)
	Uz[i]=(Dss0/8/pi/E1*uy)+(Dds1*sinADE1*uz)
	
	Vy[i]=(Dss1*x[i]/8/pi/E1*vy)+(Dds0*x[i]*sinADE1*vz)	
	Vz[i]=(Dss0*x[i]/8/pi/E1*vy)+(Dds1*x[i]*sinADE1*vz)		
	
	Wy[i]=(Dss1*x[i]/8/pi/E1*wy)+(Dds0*x[i]*sinADE1*wz)
	Wz[i]=(Dss0*x[i]/8/pi/E1*wy)+(Dds1*x[i]*sinADE1*wz)

	
end

# Transform displacements from ADCS into TDCS (Do for each component)
println("Allocation not needed here")
(Vx,Wx)  =RotateObject2D(Vx,Wx,0,0,Ct,-St) #Rotate back
(Vy,Wy)  =RotateObject2D(Vy,Wy,0,0,Ct,-St) #Rotate back
(Vz,Wz)  =RotateObject2D(Vz,Wz,0,0,Ct,-St) #Rotate back

#Add in actual burgers vector
Uy=Uy.*Dss;Vy=Vy.*Dss;Wy=Wy.*Dss;
Uz=Uz.*Dds;Vz=Vz.*Dds;Wz=Wz.*Dds;

return(Ux,Uy,Uz,Vx,Vy,Vz,Wx,Wy,Wz)
end



function TransformToADCS(y,z,Dss,Dds,SideVec,TriVertex)
#Convert to dislocation coordinate system

Ct=SideVec[3];
St=SideVec[2];
P1=TriVertex[2];
P2=TriVertex[3];
# Transform coordinates of the calculation points from TDCS into ADCS
(y1,z1)  =RotateObject2D(y,z,P1,P2,Ct,St)
# Transform the in-plane slip vector components from TDCS into ADCS
(Dss1,Dds1)=RotateObject2D(Dss,Dds,0,0,Ct,St)

return(Ct,St,y1,z1,Dss1,Dds1)
end

function AngDisDisp(x,y,z,cosA,sinA,E1,E2,cosA2,sinADE1)
# AngDisDisp calculates the "incomplete" displacements (without the 
# Burgers' function contribution) associated with an angular dislocation in
# an elastic full-space.

eta = y*cosA-z*sinA;
zeta = y*sinA+z*cosA;
r = sqrt(x^2+y^2+z^2);

# Avoid complex results for the logarithmic terms
if zeta>r
	zeta=r;
end
if z>r
	z = r;
end

rMzeta=(r-zeta);
rMz=(r-z);

ux = (x*y/r/rMz-x*eta/r/rMzeta); #b8p/E1*
vx = (eta*sinA/rMzeta-y*eta/r/rMzeta+y.^2/r/rMz+E2*(cosA*log(rMzeta)-log(rMz)));#b8p/E1*
wx = (eta*cosA/rMzeta-y/r-eta*z/r/rMzeta-E2*sinA*log(rMzeta));#b8p/E1*
	
uy = (x^2*cosA/r/rMzeta-x^2/r/rMz-E2*(cosA*log(rMzeta)-log(rMz))); 	#by/8/pi/E1*
vy = (y*cosA/r/rMzeta-sinA*cosA/rMzeta-y/r/rMz);					#by*x/8/pi/E1*			
wy = (z*cosA/r/rMzeta-cosA2/rMzeta+1/r);							#by*x/8/pi/E1*
	
uz = (E2*log(rMzeta)-x.^2/r/rMzeta);								#bz*sinADE1*
vz = (sinA/rMzeta-y/r/rMzeta);										#bz*x*sinADE1*	
wz = (cosA/rMzeta-z/r/rMzeta);										#bz*x*sinADE1*


#Export individual components
return(ux[1],uy[1],uz[1],vx[1],vy[1],vz[1],wx[1],wy[1],wz[1])
end

function TDdisp_HarFunc(X,Y,Z,P1,P2,P3,Dss,Dds,Dn,nu)
# TDdisp_HarFunc calculates the harmonic function contribution to the
# displacements associated with a triangular dislocation in a half-space.
# The function cancels the surface normal tractions induced by the main and
# image dislocations.

# Calculate unit strike, dip and normal to TD vectors: 
(Vnorm,Vstrike,Vdip)=CalculateLocalTriCoords(P1,P2,P3)

println("Might not need to do this line below")
## Transform slip vector components from TDCS into EFCS
(bX,bY,bZ) = RotateObject3DNewCoords(Dn,Dss,Dds,0,0,0,Vnorm,Vstrike,Vdip);

# Calculate contribution of angular dislocation pair on each TD side 
(Ux1,Vx1,Wx1,Uy1,Vy1,Wy1,Uz1,Vz1,Wz1) = AngSetupDispFSC(X,Y,Z,bX,bY,bZ,P1,P2,nu,Vnorm,Vstrike,Vdip); # Side P1P2
(Ux2,Vx2,Wx2,Uy2,Vy2,Wy2,Uz2,Vz2,Wz2) = AngSetupDispFSC(X,Y,Z,bX,bY,bZ,P2,P3,nu,Vnorm,Vstrike,Vdip); # Side P2P3
(Ux3,Vx3,Wx3,Uy3,Vy3,Wy3,Uz3,Vz3,Wz3) = AngSetupDispFSC(X,Y,Z,bX,bY,bZ,P3,P1,nu,Vnorm,Vstrike,Vdip); # Side P3P1

println("Pass these directly into the func above and add in there.")
# Calculate total harmonic function contribution to displacements
	#Also adding in actual burgers vector
UxDn = (Ux1+Ux2+Ux3).*Dn;
UyDn = (Vx1+Vx2+Vx3).*Dn;
UzDn = (Wx1+Wx2+Wx3).*Dn;
UxDss= (Uy1+Uy2+Uy3).*Dss;
UyDss= (Vy1+Vy2+Vy3).*Dss;
UzDss= (Wy1+Wy2+Wy3).*Dss;
UxDds= (Uz1+Uz2+Uz3).*Dds;
UyDds= (Vz1+Vz2+Vz3).*Dds;
UzDds= (Wz1+Wz2+Wz3).*Dds;

return(UxDn,UyDn,UzDn,UxDss,UyDss,UzDss,UxDds,UyDds,UzDds);
end


function AngSetupDispFSC(X,Y,Z,bX,bY,bZ,PA,PB,nu,Vnorm,Vstrike,Vdip)
# AngSetupFSC calculates the Free Surface Correction to displacements 
# associated with angular dislocation pair on each TD side.

# Calculate TD side vector and the angle of the angular dislocation pair
(SideVec,eZ,beta)=CalcSideVec(PA,PB)

if abs(beta)<eps() || abs(pi-beta)<eps()
    ue = zeros(length(X),1);
    un = zeros(length(X),1);
    uv = zeros(length(X),1);
else
    (b1,b2,b3,I,y1A,y2A,y3A,y1B,y2B,y3B,ey1,ey2,ey3)=CalcSlipVectorDiscCoords(SideVec,eZ,X,Y,Z,PA,bX,bY,bZ,beta)
    
	println("Make sure there is no allocation in these")
	## Transform slip vector components from TDCS into EFCS
	(Dn1__,Dss0n_,Dds0n_) = RotateObject3DNewCoords(1.,0.,0.,0,0,0,Vnorm,Vstrike,Vdip);
	(Dn0ss,Dss1__,Dds1ss) = RotateObject3DNewCoords(0.,1.,0.,0,0,0,Vnorm,Vstrike,Vdip);
	(Dn0ds,Dss0ds,Dds1__) = RotateObject3DNewCoords(0.,0.,1.,0,0,0,Vnorm,Vstrike,Vdip);
	# Transform slip vector components from EFCS to ADCS
	(Dn1__,Dss0n_,Dds0n_)=RotateObject3DNewCoords(Dn1__,Dss0n_,Dds0n_,0,0,0,ey1,ey2,ey3)
	(Dn0ss,Dss1__,Dds1ss)=RotateObject3DNewCoords(Dn0ss,Dss1__,Dds1ss,0,0,0,ey1,ey2,ey3)
	(Dn0ds,Dss0ds,Dds1__)=RotateObject3DNewCoords(Dn0ds,Dss0ds,Dds1__,0,0,0,ey1,ey2,ey3)
		
	println("Avoid this allocation (coming from above, be wary changing this from below...)")
	#InitOutputs
	Ux  = Array{Float64}(undef, length(X),1);
	Vx  = Array{Float64}(undef, length(X),1);
	Wx  = Array{Float64}(undef, length(X),1);
	Uy  = Array{Float64}(undef, length(X),1);
	Vy  = Array{Float64}(undef, length(X),1);
	Wy  = Array{Float64}(undef, length(X),1);
	Uz  = Array{Float64}(undef, length(X),1);
	Vz  = Array{Float64}(undef, length(X),1);
	Wz  = Array{Float64}(undef, length(X),1);

	
	Iflp=.!I; #Invert the bool

	indx=findall(I);
	indxf=findall(Iflp);
	b=-pi+beta;
	sinB = sin(b);
	cosB = cos(b);
	cotB = cot(b);
	cotB2=cot(b/2);
	# Configuration I
	for i=1:length(indx)
	
		(uxA,uyA,uzA,vxA,vyA,vzA,wxA,wyA,wzA) = AngDisDispFSC(y1A[indx[i]],y2A[indx[i]],y3A[indx[i]],cosB,sinB,cotB,cotB2,b1,b2,b3,nu,-PA[3]);
		
		
		#Add mixed components together (different coords)
		Ux[indx[i]] = -((Dn1__/4/pi/(1-nu)*uxA)+(Dss0n_/4/pi/(1-nu)*uyA)+(Dds0n_/4/pi/(1-nu)*uzA))
		Uy[indx[i]] = -((Dn0ss/4/pi/(1-nu)*uxA)+(Dss1__/4/pi/(1-nu)*uyA)+(Dds1ss/4/pi/(1-nu)*uzA))
		Uz[indx[i]] = -((Dn0ds/4/pi/(1-nu)*uxA)+(Dss0ds/4/pi/(1-nu)*uyA)+(Dds1__/4/pi/(1-nu)*uzA))
		Vx[indx[i]] = -((Dn1__/4/pi/(1-nu)*vxA)+(Dss0n_/4/pi/(1-nu)*vyA)+(Dds0n_/4/pi/(1-nu)*vzA))
		Vy[indx[i]] = -((Dn0ss/4/pi/(1-nu)*vxA)+(Dss1__/4/pi/(1-nu)*vyA)+(Dds1ss/4/pi/(1-nu)*vzA))
		Vz[indx[i]] = -((Dn0ds/4/pi/(1-nu)*vxA)+(Dss0ds/4/pi/(1-nu)*vyA)+(Dds1__/4/pi/(1-nu)*vzA))
		Wx[indx[i]] = -((Dn1__/4/pi/(1-nu)*wxA)+(Dss0n_/4/pi/(1-nu)*wyA)+(Dds0n_/4/pi/(1-nu)*wzA))
		Wy[indx[i]] = -((Dn0ss/4/pi/(1-nu)*wxA)+(Dss1__/4/pi/(1-nu)*wyA)+(Dds1ss/4/pi/(1-nu)*wzA))
		Wz[indx[i]] = -((Dn0ds/4/pi/(1-nu)*wxA)+(Dss0ds/4/pi/(1-nu)*wyA)+(Dds1__/4/pi/(1-nu)*wzA))
		
		
		(uxB,uyB,uzB,vxB,vyB,vzB,wxB,wyB,wzB) = AngDisDispFSC(y1B[indx[i]],y2B[indx[i]],y3B[indx[i]],cosB,sinB,cotB,cotB2,b1,b2,b3,nu,-PB[3]);
		
		#Add mixed components together (different coords)
		Ux[indx[i]] =Ux[indx[i]] + ((Dn1__/4/pi/(1-nu)*uxB)+(Dss0n_/4/pi/(1-nu)*uyB)+(Dds0n_/4/pi/(1-nu)*uzB))
		Uy[indx[i]] =Uy[indx[i]] + ((Dn0ss/4/pi/(1-nu)*uxB)+(Dss1__/4/pi/(1-nu)*uyB)+(Dds1ss/4/pi/(1-nu)*uzB))
		Uz[indx[i]] =Uz[indx[i]] + ((Dn0ds/4/pi/(1-nu)*uxB)+(Dss0ds/4/pi/(1-nu)*uyB)+(Dds1__/4/pi/(1-nu)*uzB))
		Vx[indx[i]] =Vx[indx[i]] + ((Dn1__/4/pi/(1-nu)*vxB)+(Dss0n_/4/pi/(1-nu)*vyB)+(Dds0n_/4/pi/(1-nu)*vzB))
		Vy[indx[i]] =Vy[indx[i]] + ((Dn0ss/4/pi/(1-nu)*vxB)+(Dss1__/4/pi/(1-nu)*vyB)+(Dds1ss/4/pi/(1-nu)*vzB))
		Vz[indx[i]] =Vz[indx[i]] + ((Dn0ds/4/pi/(1-nu)*vxB)+(Dss0ds/4/pi/(1-nu)*vyB)+(Dds1__/4/pi/(1-nu)*vzB))
		Wx[indx[i]] =Wx[indx[i]] + ((Dn1__/4/pi/(1-nu)*wxB)+(Dss0n_/4/pi/(1-nu)*wyB)+(Dds0n_/4/pi/(1-nu)*wzB))
		Wy[indx[i]] =Wy[indx[i]] + ((Dn0ss/4/pi/(1-nu)*wxB)+(Dss1__/4/pi/(1-nu)*wyB)+(Dds1ss/4/pi/(1-nu)*wzB))
		Wz[indx[i]] =Wz[indx[i]] + ((Dn0ds/4/pi/(1-nu)*wxB)+(Dss0ds/4/pi/(1-nu)*wyB)+(Dds1__/4/pi/(1-nu)*wzB))
		
	end
	
	
	b=beta;
	sinB = sin(b);
	cosB = cos(b);
	cotB = cot(b);
	cotB2=cot(b/2);
	# Configuration II
	for i=1:length(indxf)
	
		(uxA,uyA,uzA,vxA,vyA,vzA,wxA,wyA,wzA) = AngDisDispFSC(y1A[indxf[i]],y2A[indxf[i]],y3A[indxf[i]],cosB,sinB,cotB,cotB2,b1,b2,b3,nu,-PA[3]);
		
		#Add mixed components together (different coords)
		Ux[indxf[i]] = -((Dn1__/4/pi/(1-nu)*uxA)+(Dss0n_/4/pi/(1-nu)*uyA)+(Dds0n_/4/pi/(1-nu)*uzA))
		Uy[indxf[i]] = -((Dn0ss/4/pi/(1-nu)*uxA)+(Dss1__/4/pi/(1-nu)*uyA)+(Dds1ss/4/pi/(1-nu)*uzA))
		Uz[indxf[i]] = -((Dn0ds/4/pi/(1-nu)*uxA)+(Dss0ds/4/pi/(1-nu)*uyA)+(Dds1__/4/pi/(1-nu)*uzA))
		Vx[indxf[i]] = -((Dn1__/4/pi/(1-nu)*vxA)+(Dss0n_/4/pi/(1-nu)*vyA)+(Dds0n_/4/pi/(1-nu)*vzA))
		Vy[indxf[i]] = -((Dn0ss/4/pi/(1-nu)*vxA)+(Dss1__/4/pi/(1-nu)*vyA)+(Dds1ss/4/pi/(1-nu)*vzA))
		Vz[indxf[i]] = -((Dn0ds/4/pi/(1-nu)*vxA)+(Dss0ds/4/pi/(1-nu)*vyA)+(Dds1__/4/pi/(1-nu)*vzA))
		Wx[indxf[i]] = -((Dn1__/4/pi/(1-nu)*wxA)+(Dss0n_/4/pi/(1-nu)*wyA)+(Dds0n_/4/pi/(1-nu)*wzA))
		Wy[indxf[i]] = -((Dn0ss/4/pi/(1-nu)*wxA)+(Dss1__/4/pi/(1-nu)*wyA)+(Dds1ss/4/pi/(1-nu)*wzA))
		Wz[indxf[i]] = -((Dn0ds/4/pi/(1-nu)*wxA)+(Dss0ds/4/pi/(1-nu)*wyA)+(Dds1__/4/pi/(1-nu)*wzA))	
		
		(uxB,uyB,uzB,vxB,vyB,vzB,wxB,wyB,wzB) = AngDisDispFSC(y1B[indxf[i]],y2B[indxf[i]],y3B[indxf[i]],cosB,sinB,cotB,cotB2,b1,b2,b3,nu,-PB[3]);
		
		#Add mixed components together (different coords)
		Ux[indxf[i]] =Ux[indxf[i]] + ((Dn1__/4/pi/(1-nu)*uxB)+(Dss0n_/4/pi/(1-nu)*uyB)+(Dds0n_/4/pi/(1-nu)*uzB))
		Uy[indxf[i]] =Uy[indxf[i]] + ((Dn0ss/4/pi/(1-nu)*uxB)+(Dss1__/4/pi/(1-nu)*uyB)+(Dds1ss/4/pi/(1-nu)*uzB))
		Uz[indxf[i]] =Uz[indxf[i]] + ((Dn0ds/4/pi/(1-nu)*uxB)+(Dss0ds/4/pi/(1-nu)*uyB)+(Dds1__/4/pi/(1-nu)*uzB))
		Vx[indxf[i]] =Vx[indxf[i]] + ((Dn1__/4/pi/(1-nu)*vxB)+(Dss0n_/4/pi/(1-nu)*vyB)+(Dds0n_/4/pi/(1-nu)*vzB))
		Vy[indxf[i]] =Vy[indxf[i]] + ((Dn0ss/4/pi/(1-nu)*vxB)+(Dss1__/4/pi/(1-nu)*vyB)+(Dds1ss/4/pi/(1-nu)*vzB))
		Vz[indxf[i]] =Vz[indxf[i]] + ((Dn0ds/4/pi/(1-nu)*vxB)+(Dss0ds/4/pi/(1-nu)*vyB)+(Dds1__/4/pi/(1-nu)*vzB))
		Wx[indxf[i]] =Wx[indxf[i]] + ((Dn1__/4/pi/(1-nu)*wxB)+(Dss0n_/4/pi/(1-nu)*wyB)+(Dds0n_/4/pi/(1-nu)*wzB))
		Wy[indxf[i]] =Wy[indxf[i]] + ((Dn0ss/4/pi/(1-nu)*wxB)+(Dss1__/4/pi/(1-nu)*wyB)+(Dds1ss/4/pi/(1-nu)*wzB))
		Wz[indxf[i]] =Wz[indxf[i]] + ((Dn0ds/4/pi/(1-nu)*wxB)+(Dss0ds/4/pi/(1-nu)*wyB)+(Dds1__/4/pi/(1-nu)*wzB))
		
    end
	

	
	#Inverse rot mat
	VxR=[ey1[1],ey2[1],ey3[1]];
	VyR=[ey1[2],ey2[2],ey3[2]];
	VzR=[ey1[3],ey2[3],ey3[3]];
	
    # Transform total Free Surface Correction to displacements from ADCS 
    # to EFCS
	println("Reduce allocation here")
	(Ux,Vx,Wx)=RotateObject3DNewCoords(Ux,Vx,Wx,0,0,0,VxR,VyR,VzR)
	(Uy,Vy,Wy)=RotateObject3DNewCoords(Uy,Vy,Wy,0,0,0,VxR,VyR,VzR)
	(Uz,Vz,Wz)=RotateObject3DNewCoords(Uz,Vz,Wz,0,0,0,VxR,VyR,VzR)


end	#if statement

return(Ux,Vx,Wx,Uy,Vy,Wy,Uz,Vz,Wz)
end


function CalcSideVec(PA,PB)
# Calculate TD side vector and the angle of the angular dislocation pair

SideVec = PB-PA;
eZ = [0;0;1];

G=-SideVec'*eZ/norm(SideVec);
beta = acos(G[1]);

return(SideVec,eZ,beta)
end

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


function AngDisDispFSC(y1,y2,y3,cosB,sinB,cotB,cotB2,b1,b2,b3,nu,a)
# AngDisDispFSC calculates the harmonic function contribution to the 
# displacements associated with an angular dislocation in an elastic 
# half-space.

#sinB = sin(beta);
#cosB = cos(beta);
#cotB = cot(beta);
#cotB2= cot(beta/2);
y3b = y3 +2*a;
z1b = y1*cosB+y3b*sinB;
z3b = -y1*sinB+y3b*cosB;
r2b = y1 ^2 +y2 ^2 +y3b^2;
rb = sqrt(r2b);

Fib = 2*atan(-y2 /(-(rb+y3b)*cotB2+y1)); # The Burgers' function

ux = (-2*(1 -nu)*(1 -2*nu)*Fib*cotB^2 +(1 -2*nu)*y2 /
    (rb+y3b)*((1 -2*nu-a/rb)*cotB-y1 /(rb+y3b)*(nu+a/rb))+(1 -2*nu)*
    y2 *cosB*cotB/(rb+z3b)*(cosB+a/rb)+a*y2 *(y3b-a)*cotB/rb^3 +y2 *
    (y3b-a)/(rb*(rb+y3b))*(-(1 -2*nu)*cotB+y1 /(rb+y3b)*(2*nu+a/rb)+
    a*y1 /rb^2)+y2 *(y3b-a)/(rb*(rb+z3b))*(cosB/(rb+z3b)*((rb*
    cosB+y3b)*((1 -2*nu)*cosB-a/rb)*cotB+2*(1 -nu)*(rb*sinB-y1)*cosB)-
    a*y3b*cosB*cotB/rb^2)); #b1/4/pi/(1 -nu)*

vx = ((1 -2*nu)*((2*(1 -nu)*cotB^2 -nu)*log(rb+y3b)-(2*
    (1 -nu)*cotB^2 +1 -2*nu)*cosB*log(rb+z3b))-(1 -2*nu)/(rb+y3b)*(y1*
    cotB*(1 -2*nu-a/rb)+nu*y3b-a+y2 ^2 /(rb+y3b)*(nu+a/rb))-(1 -2*
    nu)*z1b*cotB/(rb+z3b)*(cosB+a/rb)-a*y1 *(y3b-a)*cotB/rb^3 +
    (y3b-a)/(rb+y3b)*(-2*nu+1 /rb*((1 -2*nu)*y1*cotB-a)+y2 ^2 /(rb*
    (rb+y3b))*(2*nu+a/rb)+a*y2 ^2 /rb^3)+(y3b-a)/(rb+z3b)*(cosB^2 -
    1 /rb*((1 -2*nu)*z1b*cotB+a*cosB)+a*y3b*z1b*cotB/rb^3 -1 /(rb*
    (rb+z3b))*(y2 ^2*cosB^2 -a*z1b*cotB/rb*(rb*cosB+y3b)))); #b1/4/pi/(1 -nu)*

wx = (2*(1 -nu)*(((1 -2*nu)*Fib*cotB)+(y2 /(rb+y3b)*(2*
    nu+a/rb))-(y2*cosB/(rb+z3b)*(cosB+a/rb)))+y2 *(y3b-a)/rb*(2*
    nu/(rb+y3b)+a/rb^2)+y2 *(y3b-a)*cosB/(rb*(rb+z3b))*(1 -2*nu-
    (rb*cosB+y3b)/(rb+z3b)*(cosB+a/rb)-a*y3b/rb^2));  #b1/4/pi/(1 -nu)*

uy = ((1 -2*nu)*((2*(1 -nu)*cotB^2 +nu)*log(rb+y3b)-(2*
    (1 -nu)*cotB^2 +1)*cosB*log(rb+z3b))+(1 -2*nu)/(rb+y3b)*(-(1 -2*nu)*
    y1*cotB+nu*y3b-a+a*y1*cotB/rb+y1 ^2 /(rb+y3b)*(nu+a/rb))-(1 -2*
    nu)*cotB/(rb+z3b)*(z1b*cosB-a*(rb*sinB-y1)/(rb*cosB))-a*y1 *
    (y3b-a)*cotB/rb^3 +(y3b-a)/(rb+y3b)*(2*nu+1 /rb*((1 -2*nu)*y1*
    cotB+a)-y1 ^2 /(rb*(rb+y3b))*(2*nu+a/rb)-a*y1 ^2 /rb^3)+(y3b-a)*
    cotB/(rb+z3b)*(-cosB*sinB+a*y1 *y3b/(rb^3*cosB)+(rb*sinB-y1)/
    rb*(2*(1 -nu)*cosB-(rb*cosB+y3b)/(rb+z3b)*(1 +a/(rb*cosB))))); #b2/4/pi/(1 -nu)*
                
vy = (2*(1 -nu)*(1 -2*nu)*Fib*cotB^2 +(1 -2*nu)*y2 /
    (rb+y3b)*(-(1 -2*nu-a/rb)*cotB+y1 /(rb+y3b)*(nu+a/rb))-(1 -2*nu)*
    y2*cotB/(rb+z3b)*(1 +a/(rb*cosB))-a*y2 *(y3b-a)*cotB/rb^3 +y2 *
    (y3b-a)/(rb*(rb+y3b))*((1 -2*nu)*cotB-2*nu*y1 /(rb+y3b)-a*y1 /rb*
    (1 /rb+1 /(rb+y3b)))+y2 *(y3b-a)*cotB/(rb*(rb+z3b))*(-2*(1 -nu)*
    cosB+(rb*cosB+y3b)/(rb+z3b)*(1 +a/(rb*cosB))+a*y3b/(rb^2*cosB))); #b2/4/pi/(1 -nu)*
                
wy = (-2*(1 -nu)*(1 -2*nu)*cotB*(log(rb+y3b)-cosB*
    log(rb+z3b))-2*(1 -nu)*y1 /(rb+y3b)*(2*nu+a/rb)+2*(1 -nu)*z1b/(rb+
    z3b)*(cosB+a/rb)+(y3b-a)/rb*((1 -2*nu)*cotB-2*nu*y1 /(rb+y3b)-a*
    y1 /rb^2)-(y3b-a)/(rb+z3b)*(cosB*sinB+(rb*cosB+y3b)*cotB/rb*
    (2*(1 -nu)*cosB-(rb*cosB+y3b)/(rb+z3b))+a/rb*(sinB-y3b*z1b/
    rb^2 -z1b*(rb*cosB+y3b)/(rb*(rb+z3b)))));  #b2/4/pi/(1 -nu)*

uz = ((1 -2*nu)*(y2 /(rb+y3b)*(1 +a/rb)-y2*cosB/(rb+
    z3b)*(cosB+a/rb))-y2 *(y3b-a)/rb*(a/rb^2 +1 /(rb+y3b))+y2 *
    (y3b-a)*cosB/(rb*(rb+z3b))*((rb*cosB+y3b)/(rb+z3b)*(cosB+a/
    rb)+a*y3b/rb^2)); #b3/4/pi/(1 -nu)*
                
vz = ((1 -2*nu)*(-sinB*log(rb+z3b)-y1 /(rb+y3b)*(1 +a/
    rb)+z1b/(rb+z3b)*(cosB+a/rb))+y1 *(y3b-a)/rb*(a/rb^2 +1 /(rb+
    y3b))-(y3b-a)/(rb+z3b)*(sinB*(cosB-a/rb)+z1b/rb*(1 +a*y3b/
    rb^2)-1 /(rb*(rb+z3b))*(y2 ^2*cosB*sinB-a*z1b/rb*(rb*cosB+y3b)))); #b3/4/pi/(1 -nu)*
                
wz = (2*(1 -nu)*Fib+2*(1 -nu)*(y2*sinB/(rb+z3b)*(cosB+
    a/rb))+y2 *(y3b-a)*sinB/(rb*(rb+z3b))*(1 +(rb*cosB+y3b)/(rb+
    z3b)*(cosB+a/rb)+a*y3b/rb^2));  #b3/4/pi/(1 -nu)*
	


# u = ux[1]+uy[1]+uz[1];
# v = vx[1]+vy[1]+vz[1];
# w = wx[1]+wy[1]+wz[1];

#Export individual components
return(ux[1],uy[1],uz[1],vx[1],vy[1],vz[1],wx[1],wy[1],wz[1])

end


function TDstrainFS(X,Y,Z,P1,P2,P3,Dss,Dds,Dn,mu,lambda)
# TDstressFS 
# Calculates stresses and strains associated with a triangular dislocation 
# in an elastic full-space.

nu =lambda/(2*(mu+lambda));    #Poisson's ratio, Equation 8.28 Pollard

#Convert from global to the dislocation coordinate system
(Vnorm,Vstrike,Vdip)=CalculateLocalTriCoords(P1,P2,P3)
(p1,p2,p3,e12,e13,e23,A,B,C,casepLog,casenLog,casezLog,x,y,z)=GlobalToTDECoords(P1,P2,P3,X,Y,Z,Vnorm,Vstrike,Vdip)
xn=x[casenLog];
yn=y[casenLog];
zn=z[casenLog];
xp=x[casepLog];
yp=y[casepLog];
zp=z[casepLog];


# Calculate first angular dislocation contribution POS
(ExxDn1Tp,EyyDn1Tp,EzzDn1Tp,ExyDn1Tp,ExzDn1Tp,EyzDn1Tp,
 ExxDss1Tp,EyyDss1Tp,EzzDss1Tp,ExyDss1Tp,ExzDss1Tp,EyzDss1Tp,
 ExxDds1Tp,EyyDds1Tp,EzzDds1Tp,ExyDds1Tp,ExzDds1Tp,EyzDds1Tp) =
 TDSetupS(xp,yp,zp,A,Dn,Dss,Dds,nu,p1,-e13);
 
# Calculate second angular dislocation contribution
(ExxDn2Tp,EyyDn2Tp,EzzDn2Tp,ExyDn2Tp,ExzDn2Tp,EyzDn2Tp,
 ExxDss2Tp,EyyDss2Tp,EzzDss2Tp,ExyDss2Tp,ExzDss2Tp,EyzDss2Tp,
 ExxDds2Tp,EyyDds2Tp,EzzDds2Tp,ExyDds2Tp,ExzDds2Tp,EyzDds2Tp) =
 TDSetupS(xp,yp,zp,B,Dn,Dss,Dds,nu,p2,e12); 
 
# Calculate third angular dislocation contribution
(ExxDn3Tp,EyyDn3Tp,EzzDn3Tp,ExyDn3Tp,ExzDn3Tp,EyzDn3Tp,
 ExxDss3Tp,EyyDss3Tp,EzzDss3Tp,ExyDss3Tp,ExzDss3Tp,EyzDss3Tp,
 ExxDds3Tp,EyyDds3Tp,EzzDds3Tp,ExyDds3Tp,ExzDds3Tp,EyzDds3Tp) =
 TDSetupS(xp,yp,zp,C,Dn,Dss,Dds,nu,p3,e23);


# Calculate first angular dislocation contribution NEG
(ExxDn1Tn,EyyDn1Tn,EzzDn1Tn,ExyDn1Tn,ExzDn1Tn,EyzDn1Tn,
 ExxDss1Tn,EyyDss1Tn,EzzDss1Tn,ExyDss1Tn,ExzDss1Tn,EyzDss1Tn,
 ExxDds1Tn,EyyDds1Tn,EzzDds1Tn,ExyDds1Tn,ExzDds1Tn,EyzDds1Tn) =
 TDSetupS(xn,yn,zn,A,Dn,Dss,Dds,nu,p1,e13);
 
# Calculate second angular dislocation contribution
(ExxDn2Tn,EyyDn2Tn,EzzDn2Tn,ExyDn2Tn,ExzDn2Tn,EyzDn2Tn,
 ExxDss2Tn,EyyDss2Tn,EzzDss2Tn,ExyDss2Tn,ExzDss2Tn,EyzDss2Tn,
 ExxDds2Tn,EyyDds2Tn,EzzDds2Tn,ExyDds2Tn,ExzDds2Tn,EyzDds2Tn) = 
 TDSetupS(xn,yn,zn,B,Dn,Dss,Dds,nu,p2,-e12);
 
# Calculate third angular dislocation contribution 
(ExxDn3Tn,EyyDn3Tn,EzzDn3Tn,ExyDn3Tn,ExzDn3Tn,EyzDn3Tn,
 ExxDss3Tn,EyyDss3Tn,EzzDss3Tn,ExyDss3Tn,ExzDss3Tn,EyzDss3Tn,
 ExxDds3Tn,EyyDds3Tn,EzzDds3Tn,ExyDds3Tn,ExzDds3Tn,EyzDds3Tn) = 
 TDSetupS(xn,yn,zn,C,Dn,Dss,Dds,nu,p3,-e23);
 
println("Remove me")
# @info e13 e12 e23
# @info ExxDn1Tp[end] EyyDn1Tp[end] EzzDn1Tp[end] ExyDn1Tp[end] ExzDn1Tp[end] EyzDn1Tp[end]
# @info ExxDn2Tp[end] EyyDn2Tp[end] EzzDn2Tp[end] ExyDn2Tp[end] ExzDn2Tp[end] EyzDn2Tp[end]
# @info ExxDn3Tp[end] EyyDn3Tp[end] EzzDn3Tp[end] ExyDn3Tp[end] ExzDn3Tp[end] EyzDn3Tp[end]
# @info ExxDn1Tn[end] EyyDn1Tn[end] EzzDn1Tn[end] ExyDn1Tn[end] ExzDn1Tn[end] EyzDn1Tn[end]
# @info ExxDn2Tn[end] EyyDn2Tn[end] EzzDn2Tn[end] ExyDn2Tn[end] ExzDn2Tn[end] EyzDn2Tn[end]
# @info ExxDn3Tn[end] EyyDn3Tn[end] EzzDn3Tn[end] ExyDn3Tn[end] ExzDn3Tn[end] EyzDn3Tn[end]
# popp 
 
#Pass directly in to funcs above, no allocation
ExxDn = Array{Float64}(undef, length(x),1); 
EyyDn = Array{Float64}(undef, length(x),1); 
EzzDn = Array{Float64}(undef, length(x),1);
ExyDn = Array{Float64}(undef, length(x),1); 
ExzDn = Array{Float64}(undef, length(x),1); 
EyzDn = Array{Float64}(undef, length(x),1);
ExxDss = Array{Float64}(undef, length(x),1); 
EyyDss = Array{Float64}(undef, length(x),1); 
EzzDss = Array{Float64}(undef, length(x),1);
ExyDss = Array{Float64}(undef, length(x),1); 
ExzDss = Array{Float64}(undef, length(x),1); 
EyzDss = Array{Float64}(undef, length(x),1);
ExxDds = Array{Float64}(undef, length(x),1); 
EyyDds = Array{Float64}(undef, length(x),1); 
EzzDds = Array{Float64}(undef, length(x),1);
ExyDds = Array{Float64}(undef, length(x),1); 
ExzDds = Array{Float64}(undef, length(x),1); 
EyzDds = Array{Float64}(undef, length(x),1);
 
#Add contributions		
ExxDn[casenLog]=(ExxDn1Tn+ExxDn2Tn+ExxDn3Tn).*Dn;
EyyDn[casenLog]=(EyyDn1Tn+EyyDn2Tn+EyyDn3Tn).*Dn;
EzzDn[casenLog]=(EzzDn1Tn+EzzDn2Tn+EzzDn3Tn).*Dn;
ExyDn[casenLog]=(ExyDn1Tn+ExyDn2Tn+ExyDn3Tn).*Dn;
ExzDn[casenLog]=(ExzDn1Tn+ExzDn2Tn+ExzDn3Tn).*Dn; 
EyzDn[casenLog]=(EyzDn1Tn+EyzDn2Tn+EyzDn3Tn).*Dn;

ExxDss[casenLog]=(ExxDss1Tn+ExxDss2Tn+ExxDss3Tn).*Dss;
EyyDss[casenLog]=(EyyDss1Tn+EyyDss2Tn+EyyDss3Tn).*Dss;
EzzDss[casenLog]=(EzzDss1Tn+EzzDss2Tn+EzzDss3Tn).*Dss;
ExyDss[casenLog]=(ExyDss1Tn+ExyDss2Tn+ExyDss3Tn).*Dss;
ExzDss[casenLog]=(ExzDss1Tn+ExzDss2Tn+ExzDss3Tn).*Dss;
EyzDss[casenLog]=(EyzDss1Tn+EyzDss2Tn+EyzDss3Tn).*Dss;

ExxDds[casenLog]=(ExxDds1Tn+ExxDds2Tn+ExxDds3Tn).*Dds;
EyyDds[casenLog]=(EyyDds1Tn+EyyDds2Tn+EyyDds3Tn).*Dds;
EzzDds[casenLog]=(EzzDds1Tn+EzzDds2Tn+EzzDds3Tn).*Dds;
ExyDds[casenLog]=(ExyDds1Tn+ExyDds2Tn+ExyDds3Tn).*Dds;
ExzDds[casenLog]=(ExzDds1Tn+ExzDds2Tn+ExzDds3Tn).*Dds;
EyzDds[casenLog]=(EyzDds1Tn+EyzDds2Tn+EyzDds3Tn).*Dds;

#Add contributions		
ExxDn[casepLog]=(ExxDn1Tp+ExxDn2Tp+ExxDn3Tp).*Dn;
EyyDn[casepLog]=(EyyDn1Tp+EyyDn2Tp+EyyDn3Tp).*Dn;
EzzDn[casepLog]=(EzzDn1Tp+EzzDn2Tp+EzzDn3Tp).*Dn;
ExyDn[casepLog]=(ExyDn1Tp+ExyDn2Tp+ExyDn3Tp).*Dn;
ExzDn[casepLog]=(ExzDn1Tp+ExzDn2Tp+ExzDn3Tp).*Dn; 
EyzDn[casepLog]=(EyzDn1Tp+EyzDn2Tp+EyzDn3Tp).*Dn;

ExxDss[casepLog]=(ExxDss1Tp+ExxDss2Tp+ExxDss3Tp).*Dss;
EyyDss[casepLog]=(EyyDss1Tp+EyyDss2Tp+EyyDss3Tp).*Dss;
EzzDss[casepLog]=(EzzDss1Tp+EzzDss2Tp+EzzDss3Tp).*Dss;
ExyDss[casepLog]=(ExyDss1Tp+ExyDss2Tp+ExyDss3Tp).*Dss;
ExzDss[casepLog]=(ExzDss1Tp+ExzDss2Tp+ExzDss3Tp).*Dss;
EyzDss[casepLog]=(EyzDss1Tp+EyzDss2Tp+EyzDss3Tp).*Dss;

ExxDds[casepLog]=(ExxDds1Tp+ExxDds2Tp+ExxDds3Tp).*Dds;
EyyDds[casepLog]=(EyyDds1Tp+EyyDds2Tp+EyyDds3Tp).*Dds;
EzzDds[casepLog]=(EzzDds1Tp+EzzDds2Tp+EzzDds3Tp).*Dds;
ExyDds[casepLog]=(ExyDds1Tp+ExyDds2Tp+ExyDds3Tp).*Dds;
ExzDds[casepLog]=(ExzDds1Tp+ExzDds2Tp+ExzDds3Tp).*Dds;
EyzDds[casepLog]=(EyzDds1Tp+EyzDds2Tp+EyzDds3Tp).*Dds;
		
		
# Calculate the strain tensor components in TDCS
for i=1:length(x)
	if casezLog[i] == 1; 
		ExxDn[i] = NaN;
		EyyDn[i] = NaN;
		EzzDn[i] = NaN;
		ExyDn[i] = NaN;
		ExzDn[i] = NaN;
		EyzDn[i] = NaN;
		
		ExxDss[i] = NaN;
		EyyDss[i] = NaN;
		EzzDss[i] = NaN;
		ExyDss[i] = NaN;
		ExzDss[i] = NaN;
		EyzDss[i] = NaN;

		ExxDds[i] = NaN;
		EyyDds[i] = NaN;
		EzzDds[i] = NaN;
		ExyDds[i] = NaN;
		ExzDds[i] = NaN;
		EyzDds[i] = NaN;
	end
end





 
println("Make sure no memory reassignments here... (func below)")
# Transform the strain tensor components from TDCS into EFCS
(ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn) = TensorTransformation3D(ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn,[Vnorm Vstrike Vdip]);

(ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss) = TensorTransformation3D(ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,[Vnorm Vstrike Vdip]);

(ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds) = TensorTransformation3D(ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,[Vnorm Vstrike Vdip]);



return(	ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn,
		ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
		ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds)
 
end


function TDSetupS(x,y,z,alpha,Dn,Dss,Dds,nu,TriVertex,SideVec)
# TDSetupS transforms coordinates of the calculation points as well as 
# slip vector components from ADCS into TDCS. It then calculates the 
# strains in ADCS and transforms them into TDCS.

(Ct,St,y1,z1,Dss1,Dds0)=TransformToADCS(y,z,1.,0.,SideVec,TriVertex)
(Dss0,Dds1)=RotateObject2D(0.,1.,0,0,Ct,St)

#Init arrays
Ang=-pi+alpha;
cosA = cos(Ang);
sinA = sin(Ang);

println("Avoid this allocation")
ExxDn = Array{Float64}(undef, length(x),1);
ExxDss= Array{Float64}(undef, length(x),1);
ExxDds= Array{Float64}(undef, length(x),1);
EyyDn = Array{Float64}(undef, length(x),1);
EyyDss= Array{Float64}(undef, length(x),1);
EyyDds= Array{Float64}(undef, length(x),1);
EzzDn = Array{Float64}(undef, length(x),1);
EzzDss= Array{Float64}(undef, length(x),1);
EzzDds= Array{Float64}(undef, length(x),1);

ExyDn = Array{Float64}(undef, length(x),1);
ExyDss= Array{Float64}(undef, length(x),1);
ExyDds= Array{Float64}(undef, length(x),1);
ExzDn = Array{Float64}(undef, length(x),1);
ExzDss= Array{Float64}(undef, length(x),1);
ExzDds= Array{Float64}(undef, length(x),1);
EyzDn = Array{Float64}(undef, length(x),1);
EyzDss= Array{Float64}(undef, length(x),1);
EyzDds= Array{Float64}(undef, length(x),1);

#Extra defs out of loop to speed it up
E1=(1-nu); #Elastic cons
E2=(2*nu+1);
cosA2=cosA^2;
sinADE1=sinA/8/pi/(1-nu);

Dn1=1;
#Burger func constants
E3=Dn1/8/pi/E1;

E5Dss1=Dss1/8/pi/E1;
E5Dss0=Dss0/8/pi/E1;




# Calculate strains associated with an angular dislocation in ADCS
for i=1:length(x)

	(exxDn,exxDss,exxDds,
	 eyyDn,eyyDss,eyyDds,
	 ezzDn,ezzDss,ezzDds,
	 exyDn,exyDss,exyDds,
	 exzDn,exzDss,exzDds,
	 eyzDn,eyzDss,eyzDds,
	 rFi_rx,rFi_ry,rFi_rz) =
	 AngDisStrain(x[i],y1[i],z1[i],cosA,sinA,Dn,Dss1[1],Dds1[1],nu,E1,E2,cosA2,sinADE1)

	#Opening components 
	ExxDn[i]	= Dn1*(rFi_rx)+E3* exxDn
	EyyDn[i] 	= Dn1*E3* eyyDn; 
	EzzDn[i] 	= Dn1*E3* ezzDn;
	ExyDn[i] 	= Dn1*(rFi_ry)/2-E3* exyDn; 
	ExzDn[i] 	= Dn1*(rFi_rz)/2-E3* exzDn;
	EyzDn[i]  	= Dn1*E3* eyzDn;
	
	#Some constants
	E4Dss1=Dss1*x[i]/8/pi/E1;
	E4Dss0=Dss0*x[i]/8/pi/E1;
	
	#Comp mixed components (in the current coords)
	ExxDss[i] 	= (-E4Dss1* exxDss)  +  (Dds0*x[i]*sinADE1* exxDds);
	ExxDds[i] 	= (-E4Dss0* exxDss)  +  (Dds1*x[i]*sinADE1* exxDds);
	
	EyyDss[i] 	= (Dss1*(rFi_ry)-E4Dss1* eyyDss)  +  (Dds0*x[i]*sinADE1* eyyDds);
	EyyDds[i] 	= (Dss0*(rFi_ry)-E4Dss0* eyyDss)  +  (Dds1*x[i]*sinADE1* eyyDds);
	
	EzzDss[i] 	= (-E4Dss1* ezzDss) + (Dds0*(rFi_rz)+Dds0*x[i]*sinADE1* ezzDds); 
	EzzDds[i] 	= (-E4Dss0* ezzDss) + (Dds1*(rFi_rz)+Dds1*x[i]*sinADE1* ezzDds); 
	
	ExyDss[i] 	= (Dss1*(rFi_rx)/2+E5Dss1* exyDss)  -  (Dds0*sinADE1* exyDds);	
	ExyDds[i] 	= (Dss0*(rFi_rx)/2+E5Dss0* exyDss)  -  (Dds1*sinADE1* exyDds);	
	
	ExzDss[i] 	= (E5Dss1* exzDss) + (Dds0*(rFi_rx)/2-Dds0*sinADE1*exzDds) ; 
	ExzDds[i] 	= (E5Dss0* exzDss) + (Dds1*(rFi_rx)/2-Dds1*sinADE1*exzDds);

	EyzDss[i] 	= (Dss1*(rFi_rz)/2-E4Dss1* eyzDss) + (Dds0*(rFi_ry)/2-Dds0*x[i]*sinADE1* eyzDds);
	EyzDds[i] 	= (Dss0*(rFi_rz)/2-E4Dss0* eyzDss) + (Dds1*(rFi_ry)/2-Dds1*x[i]*sinADE1* eyzDds);		  

end	

# @info x[1] y1[1] z1[1] cosA sinA Dn Dss1 Dds1 Dss0 Dds0 nu E1 E2 cosA2 sinADE1 Ang
# @info ExxDn[end] EyyDn[end] EzzDn[end] ExyDn[end] ExzDn[end] EyzDn[end]
# poop

# Transform strains from ADCS into TDCS
println("Remove Allocation here")
B=[[1 0 0];[0 Ct St];[0 -St Ct]]; # 3x3 Transformation matrix

#Empty arrays for tens trans (Indiv comps)
(ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn) = 
TensorTransformation3D(ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn,B);
(ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss) = 
TensorTransformation3D(ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,B);
(ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds) = 
TensorTransformation3D(ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,B);





return(	ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn,
		ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
		ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds)
end


function AngDisStrain(x,y,z,cosA,sinA,bx,by,bz,nu,E1,E2,cosA2,sinADE1)
# AngDisStrain calculates the strains associated with an angular 
# dislocation in an elastic full-space.

eta = y*cosA-z*sinA;
zeta = y*sinA+z*cosA;


x2 = x^2;
y2 = y^2;
z2 = z^2;
r2 = x2+y2+z2;
r = sqrt(r2);
r3 = r*r2;
rmz=(r-z);
rz = r*rmz;
r2z2 = r2*rmz^2;
r3z = r3*rmz;

W = zeta-r;
W2 = W^2;
Wr = W*r;
W2r = W2*r;
Wr3 = W*r3;
W2r2 = W2*r2;

C = (r*cosA-z)/Wr;
S = (r*sinA-y)/Wr;

# Partial derivatives of the Burgers' function
rFi_rx = (eta/r/(r-zeta)-y/r/rmz)/4/pi;
rFi_ry = (x/r/rmz-cosA*x/r/(r-zeta))/4/pi;
rFi_rz = (sinA*x/r/(r-zeta))/4/pi;

#Split up parts, commented are Mehdis original parts. 

#Exx = 	
ExxBx = (eta/Wr+eta*x2/W2r2-eta*x2/Wr3+y/rz-x2*y/r2z2-x2*y/r3z); #bx*(rFi_rx)+E3*
ExxBy = ((E2/Wr+x2/W2r2-x2/Wr3)*cosA+E2/rz-x2/r2z2-x2/r3z); #-E4*
ExxBz =	(E2/Wr+x2/W2r2-x2/Wr3); #bz*x*sinADE1*
		
#Eyy = 	
EyyBx = ((1/Wr+S^2-y2/Wr3)*eta+E2*y/rz-y^3/r2z2-y^3/r3z-2*nu*cosA*S); #E3*
EyyBy = (1/rz-y2/r2z2-y2/r3z+(1/Wr+S^2-y2/Wr3)*cosA); #by*(rFi_ry)-E4*!
EyyBz = (1/Wr+S^2-y2/Wr3);#bz*x*sinADE1*

#Ezz = 
EzzBx = (eta/W/r+eta*C^2-eta*z2/Wr3+y*z/r3+2*nu*sinA*C); #E3*
EzzBy = ((1/Wr+C^2-z2/Wr3)*cosA+z/r3); #-E4*
EzzBz = (1/Wr+C^2-z2/Wr3); #bz*(rFi_rz)+bz*x*sinADE1*
	
#Exy = 	
ExyBx = (x*y2/r2z2-nu*x/rz+x*y2/r3z-nu*x*cosA/Wr+eta*x*S/Wr+eta*x*y/Wr3); #bx*(rFi_ry)/2-E3*
ExyBy = (x2*y/r2z2-nu*y/rz+x2*y/r3z+nu*cosA*S+x2*y*cosA/Wr3+x2*cosA*S/Wr); #by*(rFi_rx)/2+E5*
ExyBz =	(nu*S+x2*S/Wr+x2*y/Wr3);		#-bz*sinADE1*

#Exz = 		
ExzBx =	(-x*y/r3+nu*x*sinA/Wr+eta*x*C/Wr+eta*x*z/Wr3); #bx*(rFi_rz)/2-E3*
ExzBy = (-x2/r3+nu/r+nu*cosA*C+x2*z*cosA/Wr3+x2*cosA*C/Wr);#E5*
ExzBz = (nu*C+x2*C/Wr+x2*z/Wr3);#    bz*(rFi_rx)/2+bz*sinADE1* (INNY BIT) ;

#Eyz = 
EyzBx = (y2/r3-nu/r-nu*cosA*C+nu*sinA*S+eta*sinA*cosA/W2-eta*(y*cosA+z*sinA)/W2r+eta*y*z/W2r2-eta*y*z/Wr3);	#E3*
EyzBy = (y/r3+sinA*cosA^2/W2-cosA*(y*cosA+z*sinA)/W2r+y*z*cosA/W2r2-y*z*cosA/Wr3); #by*(rFi_rz)/2-E4*
EyzBz =	(y*z/Wr3-sinA*cosA/W2+(y*cosA+z*sinA)/W2r-y*z/W2r2); #bz*(rFi_ry)/2-bz*x*sinADE1*

return(ExxBx,ExxBy,ExxBz,
	   EyyBx,EyyBy,EyyBz,
	   EzzBx,EzzBy,EzzBz,
	   ExyBx,ExyBy,ExyBz,
	   ExzBx,ExzBy,ExzBz,
	   EyzBx,EyzBy,EyzBz,
	   rFi_rx,rFi_ry,rFi_rz)
end


function TDstrain_HarFunc(X,Y,Z,P1,P2,P3,by,bz,bx,mu,lambda,nu)
# TDstrain_HarFunc calculates the harmonic function contribution to the
# strains and stresses associated with a triangular dislocation in a 
# half-space. The function cancels the surface normal tractions induced by 
# the main and image dislocations.


# Calculate unit strike, dip and normal to TD vectors: For a horizontal TD 
# as an exception, if the normal vector points upward, the strike and dip 
# vectors point Northward and Westward, whereas if the normal vector points
# downward, the strike and dip vectors point Southward and Westward, 
# respectively.
(Vnorm,Vstrike,Vdip)=CalculateLocalTriCoords(P1,P2,P3)

## Transform slip vector components from TDCS into EFCS
(bX,bY,bZ) = RotateObject3DNewCoords(bx,by,bz,0,0,0,Vnorm,Vstrike,Vdip);

# Calculate contribution of angular dislocation pair on each TD side 
(Exx1,Eyy1,Ezz1,Exy1,Exz1,Eyz1) = AngSetupStrainFSC(X,Y,Z,bX,bY,bZ,P1,P2,mu,lambda,nu); # P1P2
(Exx2,Eyy2,Ezz2,Exy2,Exz2,Eyz2) = AngSetupStrainFSC(X,Y,Z,bX,bY,bZ,P2,P3,mu,lambda,nu); # P2P3
(Exx3,Eyy3,Ezz3,Exy3,Exz3,Eyz3) = AngSetupStrainFSC(X,Y,Z,bX,bY,bZ,P3,P1,mu,lambda,nu); # P3P1

# Calculate total harmonic function contribution to strains and stresses
Exx=Exx1.+Exx2.+Exx3;
Eyy=Eyy1.+Eyy2.+Eyy3;
Ezz=Ezz1.+Ezz2.+Ezz3;
Exy=Exy1.+Exy2.+Exy3;
Exz=Exz1.+Exz2.+Exz3;
Eyz=Eyz1.+Eyz2.+Eyz3;

return(Exx,Eyy,Ezz,Exy,Exz,Eyz);
end


function AngSetupStrainFSC(X,Y,Z,bX,bY,bZ,PA,PB,mu,lambda,nu)
# AngSetupFSC_S calculates the Free Surface Correction to strains and 
# stresses associated with angular dislocation pair on each TD side.

# Calculate TD side vector and the angle of the angular dislocation pair
(SideVec,eZ,beta)=CalcSideVec(PA,PB)

if abs(beta)<eps() || abs(pi-beta)<eps()
    Exx = zeros(length(X),1);
	Eyy = zeros(length(X),1);
	Ezz = zeros(length(X),1);
	Exy = zeros(length(X),1);
	Exz = zeros(length(X),1);
	Eyz = zeros(length(X),1);
else
    (b1,b2,b3,I,y1A,y2A,y3A,y1B,y2B,y3B,ey1,ey2,ey3)=CalcSlipVectorDiscCoords(SideVec,eZ,X,Y,Z,PA,bX,bY,bZ,beta)
	AFlip = [ey1[1] ey1[2] ey1[3] ey2[1] ey2[2]  ey2[3]  ey3[1] ey3[2]  ey3[3]]; # Transformation matrix
	
    Iflp=.!I; #Invert the bool
	indx=findall(I);
	indxf=findall(Iflp);
	
    # For singularities at surface
    v11A = zeros(length(X),1);
    v22A = zeros(length(X),1);
    v33A = zeros(length(X),1);
    v12A = zeros(length(X),1);
    v13A = zeros(length(X),1);
    v23A = zeros(length(X),1);
    
    v11B = zeros(length(X),1);
    v22B = zeros(length(X),1);
    v33B = zeros(length(X),1);
    v12B = zeros(length(X),1);
    v13B = zeros(length(X),1);
    v23B = zeros(length(X),1);
    
	#Init some vars
	b=pi-beta;
	sinB = sin(b);
	cosB = cos(b);
	cotB = cot(b);
	
    # Configuration I
	for i=1:length(indx)
		(v11A[indx[i]],v22A[indx[i]],v33A[indx[i]],v12A[indx[i]],v13A[indx[i]],v23A[indx[i]]) = AngDisStrainFSC(-y1A[indx[i]],-y2A[indx[i]],y3A[indx[i]],cosB,sinB,cotB,-b1,-b2,b3,nu,-PA[3]);
		v13A[indx[i]] = -v13A[indx[i]];
		v23A[indx[i]] = -v23A[indx[i]];
    
		(v11B[indx[i]],v22B[indx[i]],v33B[indx[i]],v12B[indx[i]],v13B[indx[i]],v23B[indx[i]]) = AngDisStrainFSC(-y1B[indx[i]],-y2B[indx[i]],y3B[indx[i]],cosB,sinB,cotB,-b1,-b2,b3,nu,-PB[3]);
		v13B[indx[i]] = -v13B[indx[i]];
		v23B[indx[i]] = -v23B[indx[i]];
	end
    
	b=beta;
	sinB = sin(b);
	cosB = cos(b);
	cotB = cot(b);
    # Configuration II
	for i=1:length(indxf)
		(v11A[indxf[i]],v22A[indxf[i]],v33A[indxf[i]],v12A[indxf[i]],v13A[indxf[i]],v23A[indxf[i]]) = AngDisStrainFSC(y1A[indxf[i]],y2A[indxf[i]],y3A[indxf[i]],cosB,sinB,cotB,b1,b2,b3,nu,-PA[3]);
		
		(v11B[indxf[i]],v22B[indxf[i]],v33B[indxf[i]],v12B[indxf[i]],v13B[indxf[i]],v23B[indxf[i]]) = AngDisStrainFSC(y1B[indxf[i]],y2B[indxf[i]],y3B[indxf[i]],cosB,sinB,cotB,b1,b2,b3,nu,-PB[3]);
	end



    # Calculate total Free Surface Correction to strains in ADCS
    v11 = v11B-v11A;
    v22 = v22B-v22A;
    v33 = v33B-v33A;
    v12 = v12B-v12A;
    v13 = v13B-v13A;
    v23 = v23B-v23A;
	

    # Transform total Free Surface Correction to strains from ADCS to EFCS
    (Exx,Eyy,Ezz,Exy,Exz,Eyz) = TensorTransformation3D(v11,v22,v33,v12,v13,v23,AFlip);
end

return(Exx,Eyy,Ezz,Exy,Exz,Eyz);
end


function AngDisStrainFSC(y1,y2,y3,cosB,sinB,cotB,b1,b2,b3,nu,a)
#AngDisStrainFSC calculates the harmonic function contribution to the
#strains associated with an angular dislocation in an elastic half-space

	
#sinB=sin(beta);
#cosB=cos(beta);
#cotB=cot(beta);
y3b=y3+2*a;
z1b=y1*cosB+y3b*sinB;
z3b=-y1*sinB+y3b*cosB;
rb2=y1^2+y2^2+y3b^2;
rb=sqrt(rb2);

W1=rb*cosB+y3b;
W2=cosB+a/rb;
W3=cosB+y3b/rb;
W4=nu+a/rb;
W5=2*nu+a/rb;
W6=rb+y3b;
W7=rb+z3b;
W8=y3+a;
W9=1+a/rb/cosB;

N1=1-2*nu;

#Partial derivatives of the Burgers' function
rFib_ry2=z1b/rb/(rb+z3b)-y1/rb/(rb+y3b);#y2=xinADCS
rFib_ry1=y2/rb/(rb+y3b)-cosB*y2/rb/(rb+z3b);#y1=yinADCS
rFib_ry3=-sinB*y2/rb/(rb+z3b);#y3=zinADCS

v11 = b1*(1/4*((-2+2*nu)*N1*rFib_ry1*cotB^2-N1*y2/W6^2*((1-W5)*cotB-
    y1/W6*W4)/rb*y1+N1*y2/W6*(a/rb^3*y1*cotB-1/W6*W4+y1^2/
    W6^2*W4/rb+y1^2/W6*a/rb^3)-N1*y2*cosB*cotB/W7^2*W2*(y1/
    rb-sinB)-N1*y2*cosB*cotB/W7*a/rb^3*y1-3*a*y2*W8*cotB/rb^5*
    y1-y2*W8/rb^3/W6*(-N1*cotB+y1/W6*W5+a*y1/rb2)*y1-y2*W8/
    rb2/W6^2*(-N1*cotB+y1/W6*W5+a*y1/rb2)*y1+y2*W8/rb/W6*
    (1/W6*W5-y1^2/W6^2*W5/rb-y1^2/W6*a/rb^3+a/rb2-2*a*y1^
    2/rb2^2)-y2*W8/rb^3/W7*(cosB/W7*(W1*(N1*cosB-a/rb)*cotB+
    (2-2*nu)*(rb*sinB-y1)*cosB)-a*y3b*cosB*cotB/rb2)*y1-y2*W8/rb/
    W7^2*(cosB/W7*(W1*(N1*cosB-a/rb)*cotB+(2-2*nu)*(rb*sinB-y1)*
    cosB)-a*y3b*cosB*cotB/rb2)*(y1/rb-sinB)+y2*W8/rb/W7*(-cosB/
    W7^2*(W1*(N1*cosB-a/rb)*cotB+(2-2*nu)*(rb*sinB-y1)*cosB)*(y1/
    rb-sinB)+cosB/W7*(1/rb*cosB*y1*(N1*cosB-a/rb)*cotB+W1*a/rb^
    3*y1*cotB+(2-2*nu)*(1/rb*sinB*y1-1)*cosB)+2*a*y3b*cosB*cotB/
    rb2^2*y1))/pi/(1-nu))+
    b2*(1/4*(N1*(((2-2*nu)*cotB^2+nu)/rb*y1/W6-((2-2*nu)*cotB^2+1)*
    cosB*(y1/rb-sinB)/W7)-N1/W6^2*(-N1*y1*cotB+nu*y3b-a+a*y1*
    cotB/rb+y1^2/W6*W4)/rb*y1+N1/W6*(-N1*cotB+a*cotB/rb-a*
    y1^2*cotB/rb^3+2*y1/W6*W4-y1^3/W6^2*W4/rb-y1^3/W6*a/
    rb^3)+N1*cotB/W7^2*(z1b*cosB-a*(rb*sinB-y1)/rb/cosB)*(y1/
    rb-sinB)-N1*cotB/W7*(cosB^2-a*(1/rb*sinB*y1-1)/rb/cosB+a*
    (rb*sinB-y1)/rb^3/cosB*y1)-a*W8*cotB/rb^3+3*a*y1^2*W8*
    cotB/rb^5-W8/W6^2*(2*nu+1/rb*(N1*y1*cotB+a)-y1^2/rb/W6*
    W5-a*y1^2/rb^3)/rb*y1+W8/W6*(-1/rb^3*(N1*y1*cotB+a)*y1+
    1/rb*N1*cotB-2*y1/rb/W6*W5+y1^3/rb^3/W6*W5+y1^3/rb2/
    W6^2*W5+y1^3/rb2^2/W6*a-2*a/rb^3*y1+3*a*y1^3/rb^5)-W8*
    cotB/W7^2*(-cosB*sinB+a*y1*y3b/rb^3/cosB+(rb*sinB-y1)/rb*
    ((2-2*nu)*cosB-W1/W7*W9))*(y1/rb-sinB)+W8*cotB/W7*(a*y3b/
    rb^3/cosB-3*a*y1^2*y3b/rb^5/cosB+(1/rb*sinB*y1-1)/rb*
    ((2-2*nu)*cosB-W1/W7*W9)-(rb*sinB-y1)/rb^3*((2-2*nu)*cosB-W1/
    W7*W9)*y1+(rb*sinB-y1)/rb*(-1/rb*cosB*y1/W7*W9+W1/W7^2*
    W9*(y1/rb-sinB)+W1/W7*a/rb^3/cosB*y1)))/pi/(1-nu))+
    b3*(1/4*(N1*(-y2/W6^2*(1+a/rb)/rb*y1-y2/W6*a/rb^3*y1+y2*
    cosB/W7^2*W2*(y1/rb-sinB)+y2*cosB/W7*a/rb^3*y1)+y2*W8/
    rb^3*(a/rb2+1/W6)*y1-y2*W8/rb*(-2*a/rb2^2*y1-1/W6^2/
    rb*y1)-y2*W8*cosB/rb^3/W7*(W1/W7*W2+a*y3b/rb2)*y1-y2*W8*
    cosB/rb/W7^2*(W1/W7*W2+a*y3b/rb2)*(y1/rb-sinB)+y2*W8*
    cosB/rb/W7*(1/rb*cosB*y1/W7*W2-W1/W7^2*W2*(y1/rb-sinB)-
    W1/W7*a/rb^3*y1-2*a*y3b/rb2^2*y1))/pi/(1-nu));
	
v22 = b1*(1/4*(N1*(((2-2*nu)*cotB^2-nu)/rb*y2/W6-((2-2*nu)*cotB^2+1-
    2*nu)*cosB/rb*y2/W7)+N1/W6^2*(y1*cotB*(1-W5)+nu*y3b-a+y2^
    2/W6*W4)/rb*y2-N1/W6*(a*y1*cotB/rb^3*y2+2*y2/W6*W4-y2^
    3/W6^2*W4/rb-y2^3/W6*a/rb^3)+N1*z1b*cotB/W7^2*W2/rb*
    y2+N1*z1b*cotB/W7*a/rb^3*y2+3*a*y2*W8*cotB/rb^5*y1-W8/
    W6^2*(-2*nu+1/rb*(N1*y1*cotB-a)+y2^2/rb/W6*W5+a*y2^2/
    rb^3)/rb*y2+W8/W6*(-1/rb^3*(N1*y1*cotB-a)*y2+2*y2/rb/
    W6*W5-y2^3/rb^3/W6*W5-y2^3/rb2/W6^2*W5-y2^3/rb2^2/W6*
    a+2*a/rb^3*y2-3*a*y2^3/rb^5)-W8/W7^2*(cosB^2-1/rb*(N1*
    z1b*cotB+a*cosB)+a*y3b*z1b*cotB/rb^3-1/rb/W7*(y2^2*cosB^2-
    a*z1b*cotB/rb*W1))/rb*y2+W8/W7*(1/rb^3*(N1*z1b*cotB+a*
    cosB)*y2-3*a*y3b*z1b*cotB/rb^5*y2+1/rb^3/W7*(y2^2*cosB^2-
    a*z1b*cotB/rb*W1)*y2+1/rb2/W7^2*(y2^2*cosB^2-a*z1b*cotB/
    rb*W1)*y2-1/rb/W7*(2*y2*cosB^2+a*z1b*cotB/rb^3*W1*y2-a*
    z1b*cotB/rb2*cosB*y2)))/pi/(1-nu))+
    b2*(1/4*((2-2*nu)*N1*rFib_ry2*cotB^2+N1/W6*((W5-1)*cotB+y1/W6*
    W4)-N1*y2^2/W6^2*((W5-1)*cotB+y1/W6*W4)/rb+N1*y2/W6*(-a/
    rb^3*y2*cotB-y1/W6^2*W4/rb*y2-y2/W6*a/rb^3*y1)-N1*cotB/
    W7*W9+N1*y2^2*cotB/W7^2*W9/rb+N1*y2^2*cotB/W7*a/rb^3/
    cosB-a*W8*cotB/rb^3+3*a*y2^2*W8*cotB/rb^5+W8/rb/W6*(N1*
    cotB-2*nu*y1/W6-a*y1/rb*(1/rb+1/W6))-y2^2*W8/rb^3/W6*
    (N1*cotB-2*nu*y1/W6-a*y1/rb*(1/rb+1/W6))-y2^2*W8/rb2/W6^
    2*(N1*cotB-2*nu*y1/W6-a*y1/rb*(1/rb+1/W6))+y2*W8/rb/W6*
    (2*nu*y1/W6^2/rb*y2+a*y1/rb^3*(1/rb+1/W6)*y2-a*y1/rb*
    (-1/rb^3*y2-1/W6^2/rb*y2))+W8*cotB/rb/W7*((-2+2*nu)*cosB+
    W1/W7*W9+a*y3b/rb2/cosB)-y2^2*W8*cotB/rb^3/W7*((-2+2*nu)*
    cosB+W1/W7*W9+a*y3b/rb2/cosB)-y2^2*W8*cotB/rb2/W7^2*((-2+
    2*nu)*cosB+W1/W7*W9+a*y3b/rb2/cosB)+y2*W8*cotB/rb/W7*(1/
    rb*cosB*y2/W7*W9-W1/W7^2*W9/rb*y2-W1/W7*a/rb^3/cosB*y2-
    2*a*y3b/rb2^2/cosB*y2))/pi/(1-nu))+
    b3*(1/4*(N1*(-sinB/rb*y2/W7+y2/W6^2*(1+a/rb)/rb*y1+y2/W6*
    a/rb^3*y1-z1b/W7^2*W2/rb*y2-z1b/W7*a/rb^3*y2)-y2*W8/
    rb^3*(a/rb2+1/W6)*y1+y1*W8/rb*(-2*a/rb2^2*y2-1/W6^2/
    rb*y2)+W8/W7^2*(sinB*(cosB-a/rb)+z1b/rb*(1+a*y3b/rb2)-1/
    rb/W7*(y2^2*cosB*sinB-a*z1b/rb*W1))/rb*y2-W8/W7*(sinB*a/
    rb^3*y2-z1b/rb^3*(1+a*y3b/rb2)*y2-2*z1b/rb^5*a*y3b*y2+
    1/rb^3/W7*(y2^2*cosB*sinB-a*z1b/rb*W1)*y2+1/rb2/W7^2*
    (y2^2*cosB*sinB-a*z1b/rb*W1)*y2-1/rb/W7*(2*y2*cosB*sinB+a*
    z1b/rb^3*W1*y2-a*z1b/rb2*cosB*y2)))/pi/(1-nu));

v33 = b1*(1/4*((2-2*nu)*(N1*rFib_ry3*cotB-y2/W6^2*W5*(y3b/rb+1)-
    1/2*y2/W6*a/rb^3*2*y3b+y2*cosB/W7^2*W2*W3+1/2*y2*cosB/W7*
    a/rb^3*2*y3b)+y2/rb*(2*nu/W6+a/rb2)-1/2*y2*W8/rb^3*(2*
    nu/W6+a/rb2)*2*y3b+y2*W8/rb*(-2*nu/W6^2*(y3b/rb+1)-a/
    rb2^2*2*y3b)+y2*cosB/rb/W7*(1-2*nu-W1/W7*W2-a*y3b/rb2)-
    1/2*y2*W8*cosB/rb^3/W7*(1-2*nu-W1/W7*W2-a*y3b/rb2)*2*
    y3b-y2*W8*cosB/rb/W7^2*(1-2*nu-W1/W7*W2-a*y3b/rb2)*W3+y2*
    W8*cosB/rb/W7*(-(cosB*y3b/rb+1)/W7*W2+W1/W7^2*W2*W3+1/2*
    W1/W7*a/rb^3*2*y3b-a/rb2+a*y3b/rb2^2*2*y3b))/pi/(1-nu))+
    b2*(1/4*((-2+2*nu)*N1*cotB*((y3b/rb+1)/W6-cosB*W3/W7)+(2-2*nu)*
    y1/W6^2*W5*(y3b/rb+1)+1/2*(2-2*nu)*y1/W6*a/rb^3*2*y3b+(2-
    2*nu)*sinB/W7*W2-(2-2*nu)*z1b/W7^2*W2*W3-1/2*(2-2*nu)*z1b/
    W7*a/rb^3*2*y3b+1/rb*(N1*cotB-2*nu*y1/W6-a*y1/rb2)-1/2*
    W8/rb^3*(N1*cotB-2*nu*y1/W6-a*y1/rb2)*2*y3b+W8/rb*(2*nu*
    y1/W6^2*(y3b/rb+1)+a*y1/rb2^2*2*y3b)-1/W7*(cosB*sinB+W1*
    cotB/rb*((2-2*nu)*cosB-W1/W7)+a/rb*(sinB-y3b*z1b/rb2-z1b*
    W1/rb/W7))+W8/W7^2*(cosB*sinB+W1*cotB/rb*((2-2*nu)*cosB-W1/
    W7)+a/rb*(sinB-y3b*z1b/rb2-z1b*W1/rb/W7))*W3-W8/W7*((cosB*
    y3b/rb+1)*cotB/rb*((2-2*nu)*cosB-W1/W7)-1/2*W1*cotB/rb^3*
    ((2-2*nu)*cosB-W1/W7)*2*y3b+W1*cotB/rb*(-(cosB*y3b/rb+1)/W7+
    W1/W7^2*W3)-1/2*a/rb^3*(sinB-y3b*z1b/rb2-z1b*W1/rb/W7)*
    2*y3b+a/rb*(-z1b/rb2-y3b*sinB/rb2+y3b*z1b/rb2^2*2*y3b-
    sinB*W1/rb/W7-z1b*(cosB*y3b/rb+1)/rb/W7+1/2*z1b*W1/rb^3/
    W7*2*y3b+z1b*W1/rb/W7^2*W3)))/pi/(1-nu))+
    b3*(1/4*((2-2*nu)*rFib_ry3-(2-2*nu)*y2*sinB/W7^2*W2*W3-1/2*
    (2-2*nu)*y2*sinB/W7*a/rb^3*2*y3b+y2*sinB/rb/W7*(1+W1/W7*
    W2+a*y3b/rb2)-1/2*y2*W8*sinB/rb^3/W7*(1+W1/W7*W2+a*y3b/
    rb2)*2*y3b-y2*W8*sinB/rb/W7^2*(1+W1/W7*W2+a*y3b/rb2)*W3+
    y2*W8*sinB/rb/W7*((cosB*y3b/rb+1)/W7*W2-W1/W7^2*W2*W3-
    1/2*W1/W7*a/rb^3*2*y3b+a/rb2-a*y3b/rb2^2*2*y3b))/pi/(1-nu));

v12 = b1/2*(1/4*((-2+2*nu)*N1*rFib_ry2*cotB^2+N1/W6*((1-W5)*cotB-y1/
    W6*W4)-N1*y2^2/W6^2*((1-W5)*cotB-y1/W6*W4)/rb+N1*y2/W6*
    (a/rb^3*y2*cotB+y1/W6^2*W4/rb*y2+y2/W6*a/rb^3*y1)+N1*
    cosB*cotB/W7*W2-N1*y2^2*cosB*cotB/W7^2*W2/rb-N1*y2^2*cosB*
    cotB/W7*a/rb^3+a*W8*cotB/rb^3-3*a*y2^2*W8*cotB/rb^5+W8/
    rb/W6*(-N1*cotB+y1/W6*W5+a*y1/rb2)-y2^2*W8/rb^3/W6*(-N1*
    cotB+y1/W6*W5+a*y1/rb2)-y2^2*W8/rb2/W6^2*(-N1*cotB+y1/
    W6*W5+a*y1/rb2)+y2*W8/rb/W6*(-y1/W6^2*W5/rb*y2-y2/W6*
    a/rb^3*y1-2*a*y1/rb2^2*y2)+W8/rb/W7*(cosB/W7*(W1*(N1*
    cosB-a/rb)*cotB+(2-2*nu)*(rb*sinB-y1)*cosB)-a*y3b*cosB*cotB/
    rb2)-y2^2*W8/rb^3/W7*(cosB/W7*(W1*(N1*cosB-a/rb)*cotB+(2-
    2*nu)*(rb*sinB-y1)*cosB)-a*y3b*cosB*cotB/rb2)-y2^2*W8/rb2/
    W7^2*(cosB/W7*(W1*(N1*cosB-a/rb)*cotB+(2-2*nu)*(rb*sinB-y1)*
    cosB)-a*y3b*cosB*cotB/rb2)+y2*W8/rb/W7*(-cosB/W7^2*(W1*
    (N1*cosB-a/rb)*cotB+(2-2*nu)*(rb*sinB-y1)*cosB)/rb*y2+cosB/
    W7*(1/rb*cosB*y2*(N1*cosB-a/rb)*cotB+W1*a/rb^3*y2*cotB+(2-2*
    nu)/rb*sinB*y2*cosB)+2*a*y3b*cosB*cotB/rb2^2*y2))/pi/(1-nu))+
    b2/2*(1/4*(N1*(((2-2*nu)*cotB^2+nu)/rb*y2/W6-((2-2*nu)*cotB^2+1)*
    cosB/rb*y2/W7)-N1/W6^2*(-N1*y1*cotB+nu*y3b-a+a*y1*cotB/rb+
    y1^2/W6*W4)/rb*y2+N1/W6*(-a*y1*cotB/rb^3*y2-y1^2/W6^
    2*W4/rb*y2-y1^2/W6*a/rb^3*y2)+N1*cotB/W7^2*(z1b*cosB-a*
    (rb*sinB-y1)/rb/cosB)/rb*y2-N1*cotB/W7*(-a/rb2*sinB*y2/
    cosB+a*(rb*sinB-y1)/rb^3/cosB*y2)+3*a*y2*W8*cotB/rb^5*y1-
    W8/W6^2*(2*nu+1/rb*(N1*y1*cotB+a)-y1^2/rb/W6*W5-a*y1^2/
    rb^3)/rb*y2+W8/W6*(-1/rb^3*(N1*y1*cotB+a)*y2+y1^2/rb^
    3/W6*W5*y2+y1^2/rb2/W6^2*W5*y2+y1^2/rb2^2/W6*a*y2+3*
    a*y1^2/rb^5*y2)-W8*cotB/W7^2*(-cosB*sinB+a*y1*y3b/rb^3/
    cosB+(rb*sinB-y1)/rb*((2-2*nu)*cosB-W1/W7*W9))/rb*y2+W8*cotB/
    W7*(-3*a*y1*y3b/rb^5/cosB*y2+1/rb2*sinB*y2*((2-2*nu)*cosB-
    W1/W7*W9)-(rb*sinB-y1)/rb^3*((2-2*nu)*cosB-W1/W7*W9)*y2+(rb*
    sinB-y1)/rb*(-1/rb*cosB*y2/W7*W9+W1/W7^2*W9/rb*y2+W1/W7*
    a/rb^3/cosB*y2)))/pi/(1-nu))+
    b3/2*(1/4*(N1*(1/W6*(1+a/rb)-y2^2/W6^2*(1+a/rb)/rb-y2^2/
    W6*a/rb^3-cosB/W7*W2+y2^2*cosB/W7^2*W2/rb+y2^2*cosB/W7*
    a/rb^3)-W8/rb*(a/rb2+1/W6)+y2^2*W8/rb^3*(a/rb2+1/W6)-
    y2*W8/rb*(-2*a/rb2^2*y2-1/W6^2/rb*y2)+W8*cosB/rb/W7*
    (W1/W7*W2+a*y3b/rb2)-y2^2*W8*cosB/rb^3/W7*(W1/W7*W2+a*
    y3b/rb2)-y2^2*W8*cosB/rb2/W7^2*(W1/W7*W2+a*y3b/rb2)+y2*
    W8*cosB/rb/W7*(1/rb*cosB*y2/W7*W2-W1/W7^2*W2/rb*y2-W1/
    W7*a/rb^3*y2-2*a*y3b/rb2^2*y2))/pi/(1-nu))+
    b1/2*(1/4*(N1*(((2-2*nu)*cotB^2-nu)/rb*y1/W6-((2-2*nu)*cotB^2+1-
    2*nu)*cosB*(y1/rb-sinB)/W7)+N1/W6^2*(y1*cotB*(1-W5)+nu*y3b-
    a+y2^2/W6*W4)/rb*y1-N1/W6*((1-W5)*cotB+a*y1^2*cotB/rb^3-
    y2^2/W6^2*W4/rb*y1-y2^2/W6*a/rb^3*y1)-N1*cosB*cotB/W7*
    W2+N1*z1b*cotB/W7^2*W2*(y1/rb-sinB)+N1*z1b*cotB/W7*a/rb^
    3*y1-a*W8*cotB/rb^3+3*a*y1^2*W8*cotB/rb^5-W8/W6^2*(-2*
    nu+1/rb*(N1*y1*cotB-a)+y2^2/rb/W6*W5+a*y2^2/rb^3)/rb*
    y1+W8/W6*(-1/rb^3*(N1*y1*cotB-a)*y1+1/rb*N1*cotB-y2^2/
    rb^3/W6*W5*y1-y2^2/rb2/W6^2*W5*y1-y2^2/rb2^2/W6*a*y1-
    3*a*y2^2/rb^5*y1)-W8/W7^2*(cosB^2-1/rb*(N1*z1b*cotB+a*
    cosB)+a*y3b*z1b*cotB/rb^3-1/rb/W7*(y2^2*cosB^2-a*z1b*cotB/
    rb*W1))*(y1/rb-sinB)+W8/W7*(1/rb^3*(N1*z1b*cotB+a*cosB)*
    y1-1/rb*N1*cosB*cotB+a*y3b*cosB*cotB/rb^3-3*a*y3b*z1b*cotB/
    rb^5*y1+1/rb^3/W7*(y2^2*cosB^2-a*z1b*cotB/rb*W1)*y1+1/
    rb/W7^2*(y2^2*cosB^2-a*z1b*cotB/rb*W1)*(y1/rb-sinB)-1/rb/
    W7*(-a*cosB*cotB/rb*W1+a*z1b*cotB/rb^3*W1*y1-a*z1b*cotB/
    rb2*cosB*y1)))/pi/(1-nu))+
    b2/2*(1/4*((2-2*nu)*N1*rFib_ry1*cotB^2-N1*y2/W6^2*((W5-1)*cotB+
    y1/W6*W4)/rb*y1+N1*y2/W6*(-a/rb^3*y1*cotB+1/W6*W4-y1^
    2/W6^2*W4/rb-y1^2/W6*a/rb^3)+N1*y2*cotB/W7^2*W9*(y1/
    rb-sinB)+N1*y2*cotB/W7*a/rb^3/cosB*y1+3*a*y2*W8*cotB/rb^
    5*y1-y2*W8/rb^3/W6*(N1*cotB-2*nu*y1/W6-a*y1/rb*(1/rb+1/
    W6))*y1-y2*W8/rb2/W6^2*(N1*cotB-2*nu*y1/W6-a*y1/rb*(1/
    rb+1/W6))*y1+y2*W8/rb/W6*(-2*nu/W6+2*nu*y1^2/W6^2/rb-a/
    rb*(1/rb+1/W6)+a*y1^2/rb^3*(1/rb+1/W6)-a*y1/rb*(-1/
    rb^3*y1-1/W6^2/rb*y1))-y2*W8*cotB/rb^3/W7*((-2+2*nu)*
    cosB+W1/W7*W9+a*y3b/rb2/cosB)*y1-y2*W8*cotB/rb/W7^2*((-2+
    2*nu)*cosB+W1/W7*W9+a*y3b/rb2/cosB)*(y1/rb-sinB)+y2*W8*
    cotB/rb/W7*(1/rb*cosB*y1/W7*W9-W1/W7^2*W9*(y1/rb-sinB)-
    W1/W7*a/rb^3/cosB*y1-2*a*y3b/rb2^2/cosB*y1))/pi/(1-nu))+
    b3/2*(1/4*(N1*(-sinB*(y1/rb-sinB)/W7-1/W6*(1+a/rb)+y1^2/W6^
    2*(1+a/rb)/rb+y1^2/W6*a/rb^3+cosB/W7*W2-z1b/W7^2*W2*
    (y1/rb-sinB)-z1b/W7*a/rb^3*y1)+W8/rb*(a/rb2+1/W6)-y1^2*
    W8/rb^3*(a/rb2+1/W6)+y1*W8/rb*(-2*a/rb2^2*y1-1/W6^2/
    rb*y1)+W8/W7^2*(sinB*(cosB-a/rb)+z1b/rb*(1+a*y3b/rb2)-1/
    rb/W7*(y2^2*cosB*sinB-a*z1b/rb*W1))*(y1/rb-sinB)-W8/W7*
    (sinB*a/rb^3*y1+cosB/rb*(1+a*y3b/rb2)-z1b/rb^3*(1+a*y3b/
    rb2)*y1-2*z1b/rb^5*a*y3b*y1+1/rb^3/W7*(y2^2*cosB*sinB-a*
    z1b/rb*W1)*y1+1/rb/W7^2*(y2^2*cosB*sinB-a*z1b/rb*W1)*
    (y1/rb-sinB)-1/rb/W7*(-a*cosB/rb*W1+a*z1b/rb^3*W1*y1-a*
    z1b/rb2*cosB*y1)))/pi/(1-nu));

v13 = b1/2*(1/4*((-2+2*nu)*N1*rFib_ry3*cotB^2-N1*y2/W6^2*((1-W5)*
    cotB-y1/W6*W4)*(y3b/rb+1)+N1*y2/W6*(1/2*a/rb^3*2*y3b*cotB+
    y1/W6^2*W4*(y3b/rb+1)+1/2*y1/W6*a/rb^3*2*y3b)-N1*y2*cosB*
    cotB/W7^2*W2*W3-1/2*N1*y2*cosB*cotB/W7*a/rb^3*2*y3b+a/
    rb^3*y2*cotB-3/2*a*y2*W8*cotB/rb^5*2*y3b+y2/rb/W6*(-N1*
    cotB+y1/W6*W5+a*y1/rb2)-1/2*y2*W8/rb^3/W6*(-N1*cotB+y1/
    W6*W5+a*y1/rb2)*2*y3b-y2*W8/rb/W6^2*(-N1*cotB+y1/W6*W5+
    a*y1/rb2)*(y3b/rb+1)+y2*W8/rb/W6*(-y1/W6^2*W5*(y3b/rb+
    1)-1/2*y1/W6*a/rb^3*2*y3b-a*y1/rb2^2*2*y3b)+y2/rb/W7*
    (cosB/W7*(W1*(N1*cosB-a/rb)*cotB+(2-2*nu)*(rb*sinB-y1)*cosB)-
    a*y3b*cosB*cotB/rb2)-1/2*y2*W8/rb^3/W7*(cosB/W7*(W1*(N1*
    cosB-a/rb)*cotB+(2-2*nu)*(rb*sinB-y1)*cosB)-a*y3b*cosB*cotB/
    rb2)*2*y3b-y2*W8/rb/W7^2*(cosB/W7*(W1*(N1*cosB-a/rb)*cotB+
    (2-2*nu)*(rb*sinB-y1)*cosB)-a*y3b*cosB*cotB/rb2)*W3+y2*W8/rb/
    W7*(-cosB/W7^2*(W1*(N1*cosB-a/rb)*cotB+(2-2*nu)*(rb*sinB-y1)*
    cosB)*W3+cosB/W7*((cosB*y3b/rb+1)*(N1*cosB-a/rb)*cotB+1/2*W1*
    a/rb^3*2*y3b*cotB+1/2*(2-2*nu)/rb*sinB*2*y3b*cosB)-a*cosB*
    cotB/rb2+a*y3b*cosB*cotB/rb2^2*2*y3b))/pi/(1-nu))+
    b2/2*(1/4*(N1*(((2-2*nu)*cotB^2+nu)*(y3b/rb+1)/W6-((2-2*nu)*cotB^
    2+1)*cosB*W3/W7)-N1/W6^2*(-N1*y1*cotB+nu*y3b-a+a*y1*cotB/
    rb+y1^2/W6*W4)*(y3b/rb+1)+N1/W6*(nu-1/2*a*y1*cotB/rb^3*2*
    y3b-y1^2/W6^2*W4*(y3b/rb+1)-1/2*y1^2/W6*a/rb^3*2*y3b)+
    N1*cotB/W7^2*(z1b*cosB-a*(rb*sinB-y1)/rb/cosB)*W3-N1*cotB/
    W7*(cosB*sinB-1/2*a/rb2*sinB*2*y3b/cosB+1/2*a*(rb*sinB-y1)/
    rb^3/cosB*2*y3b)-a/rb^3*y1*cotB+3/2*a*y1*W8*cotB/rb^5*2*
    y3b+1/W6*(2*nu+1/rb*(N1*y1*cotB+a)-y1^2/rb/W6*W5-a*y1^2/
    rb^3)-W8/W6^2*(2*nu+1/rb*(N1*y1*cotB+a)-y1^2/rb/W6*W5-a*
    y1^2/rb^3)*(y3b/rb+1)+W8/W6*(-1/2/rb^3*(N1*y1*cotB+a)*2*
    y3b+1/2*y1^2/rb^3/W6*W5*2*y3b+y1^2/rb/W6^2*W5*(y3b/rb+
    1)+1/2*y1^2/rb2^2/W6*a*2*y3b+3/2*a*y1^2/rb^5*2*y3b)+
    cotB/W7*(-cosB*sinB+a*y1*y3b/rb^3/cosB+(rb*sinB-y1)/rb*((2-
    2*nu)*cosB-W1/W7*W9))-W8*cotB/W7^2*(-cosB*sinB+a*y1*y3b/rb^
    3/cosB+(rb*sinB-y1)/rb*((2-2*nu)*cosB-W1/W7*W9))*W3+W8*cotB/
    W7*(a/rb^3/cosB*y1-3/2*a*y1*y3b/rb^5/cosB*2*y3b+1/2/
    rb2*sinB*2*y3b*((2-2*nu)*cosB-W1/W7*W9)-1/2*(rb*sinB-y1)/rb^
    3*((2-2*nu)*cosB-W1/W7*W9)*2*y3b+(rb*sinB-y1)/rb*(-(cosB*y3b/
    rb+1)/W7*W9+W1/W7^2*W9*W3+1/2*W1/W7*a/rb^3/cosB*2*
    y3b)))/pi/(1-nu))+
    b3/2*(1/4*(N1*(-y2/W6^2*(1+a/rb)*(y3b/rb+1)-1/2*y2/W6*a/
    rb^3*2*y3b+y2*cosB/W7^2*W2*W3+1/2*y2*cosB/W7*a/rb^3*2*
    y3b)-y2/rb*(a/rb2+1/W6)+1/2*y2*W8/rb^3*(a/rb2+1/W6)*2*
    y3b-y2*W8/rb*(-a/rb2^2*2*y3b-1/W6^2*(y3b/rb+1))+y2*cosB/
    rb/W7*(W1/W7*W2+a*y3b/rb2)-1/2*y2*W8*cosB/rb^3/W7*(W1/
    W7*W2+a*y3b/rb2)*2*y3b-y2*W8*cosB/rb/W7^2*(W1/W7*W2+a*
    y3b/rb2)*W3+y2*W8*cosB/rb/W7*((cosB*y3b/rb+1)/W7*W2-W1/
    W7^2*W2*W3-1/2*W1/W7*a/rb^3*2*y3b+a/rb2-a*y3b/rb2^2*2*
    y3b))/pi/(1-nu))+
    b1/2*(1/4*((2-2*nu)*(N1*rFib_ry1*cotB-y1/W6^2*W5/rb*y2-y2/W6*
    a/rb^3*y1+y2*cosB/W7^2*W2*(y1/rb-sinB)+y2*cosB/W7*a/rb^
    3*y1)-y2*W8/rb^3*(2*nu/W6+a/rb2)*y1+y2*W8/rb*(-2*nu/W6^
    2/rb*y1-2*a/rb2^2*y1)-y2*W8*cosB/rb^3/W7*(1-2*nu-W1/W7*
    W2-a*y3b/rb2)*y1-y2*W8*cosB/rb/W7^2*(1-2*nu-W1/W7*W2-a*
    y3b/rb2)*(y1/rb-sinB)+y2*W8*cosB/rb/W7*(-1/rb*cosB*y1/W7*
    W2+W1/W7^2*W2*(y1/rb-sinB)+W1/W7*a/rb^3*y1+2*a*y3b/rb2^
    2*y1))/pi/(1-nu))+
    b2/2*(1/4*((-2+2*nu)*N1*cotB*(1/rb*y1/W6-cosB*(y1/rb-sinB)/W7)-
    (2-2*nu)/W6*W5+(2-2*nu)*y1^2/W6^2*W5/rb+(2-2*nu)*y1^2/W6*
    a/rb^3+(2-2*nu)*cosB/W7*W2-(2-2*nu)*z1b/W7^2*W2*(y1/rb-
    sinB)-(2-2*nu)*z1b/W7*a/rb^3*y1-W8/rb^3*(N1*cotB-2*nu*y1/
    W6-a*y1/rb2)*y1+W8/rb*(-2*nu/W6+2*nu*y1^2/W6^2/rb-a/rb2+
    2*a*y1^2/rb2^2)+W8/W7^2*(cosB*sinB+W1*cotB/rb*((2-2*nu)*
    cosB-W1/W7)+a/rb*(sinB-y3b*z1b/rb2-z1b*W1/rb/W7))*(y1/rb-
    sinB)-W8/W7*(1/rb2*cosB*y1*cotB*((2-2*nu)*cosB-W1/W7)-W1*
    cotB/rb^3*((2-2*nu)*cosB-W1/W7)*y1+W1*cotB/rb*(-1/rb*cosB*
    y1/W7+W1/W7^2*(y1/rb-sinB))-a/rb^3*(sinB-y3b*z1b/rb2-
    z1b*W1/rb/W7)*y1+a/rb*(-y3b*cosB/rb2+2*y3b*z1b/rb2^2*y1-
    cosB*W1/rb/W7-z1b/rb2*cosB*y1/W7+z1b*W1/rb^3/W7*y1+z1b*
    W1/rb/W7^2*(y1/rb-sinB))))/pi/(1-nu))+
    b3/2*(1/4*((2-2*nu)*rFib_ry1-(2-2*nu)*y2*sinB/W7^2*W2*(y1/rb-
    sinB)-(2-2*nu)*y2*sinB/W7*a/rb^3*y1-y2*W8*sinB/rb^3/W7*(1+
    W1/W7*W2+a*y3b/rb2)*y1-y2*W8*sinB/rb/W7^2*(1+W1/W7*W2+
    a*y3b/rb2)*(y1/rb-sinB)+y2*W8*sinB/rb/W7*(1/rb*cosB*y1/
    W7*W2-W1/W7^2*W2*(y1/rb-sinB)-W1/W7*a/rb^3*y1-2*a*y3b/
    rb2^2*y1))/pi/(1-nu));

v23 = b1/2*(1/4*(N1*(((2-2*nu)*cotB^2-nu)*(y3b/rb+1)/W6-((2-2*nu)*
    cotB^2+1-2*nu)*cosB*W3/W7)+N1/W6^2*(y1*cotB*(1-W5)+nu*y3b-a+
    y2^2/W6*W4)*(y3b/rb+1)-N1/W6*(1/2*a*y1*cotB/rb^3*2*y3b+
    nu-y2^2/W6^2*W4*(y3b/rb+1)-1/2*y2^2/W6*a/rb^3*2*y3b)-N1*
    sinB*cotB/W7*W2+N1*z1b*cotB/W7^2*W2*W3+1/2*N1*z1b*cotB/W7*
    a/rb^3*2*y3b-a/rb^3*y1*cotB+3/2*a*y1*W8*cotB/rb^5*2*y3b+
    1/W6*(-2*nu+1/rb*(N1*y1*cotB-a)+y2^2/rb/W6*W5+a*y2^2/
    rb^3)-W8/W6^2*(-2*nu+1/rb*(N1*y1*cotB-a)+y2^2/rb/W6*W5+
    a*y2^2/rb^3)*(y3b/rb+1)+W8/W6*(-1/2/rb^3*(N1*y1*cotB-a)*
    2*y3b-1/2*y2^2/rb^3/W6*W5*2*y3b-y2^2/rb/W6^2*W5*(y3b/
    rb+1)-1/2*y2^2/rb2^2/W6*a*2*y3b-3/2*a*y2^2/rb^5*2*y3b)+
    1/W7*(cosB^2-1/rb*(N1*z1b*cotB+a*cosB)+a*y3b*z1b*cotB/rb^
    3-1/rb/W7*(y2^2*cosB^2-a*z1b*cotB/rb*W1))-W8/W7^2*(cosB^2-
    1/rb*(N1*z1b*cotB+a*cosB)+a*y3b*z1b*cotB/rb^3-1/rb/W7*
    (y2^2*cosB^2-a*z1b*cotB/rb*W1))*W3+W8/W7*(1/2/rb^3*(N1*
    z1b*cotB+a*cosB)*2*y3b-1/rb*N1*sinB*cotB+a*z1b*cotB/rb^3+a*
    y3b*sinB*cotB/rb^3-3/2*a*y3b*z1b*cotB/rb^5*2*y3b+1/2/rb^
    3/W7*(y2^2*cosB^2-a*z1b*cotB/rb*W1)*2*y3b+1/rb/W7^2*(y2^
    2*cosB^2-a*z1b*cotB/rb*W1)*W3-1/rb/W7*(-a*sinB*cotB/rb*W1+
    1/2*a*z1b*cotB/rb^3*W1*2*y3b-a*z1b*cotB/rb*(cosB*y3b/rb+
    1))))/pi/(1-nu))+
    b2/2*(1/4*((2-2*nu)*N1*rFib_ry3*cotB^2-N1*y2/W6^2*((W5-1)*cotB+
    y1/W6*W4)*(y3b/rb+1)+N1*y2/W6*(-1/2*a/rb^3*2*y3b*cotB-y1/
    W6^2*W4*(y3b/rb+1)-1/2*y1/W6*a/rb^3*2*y3b)+N1*y2*cotB/
    W7^2*W9*W3+1/2*N1*y2*cotB/W7*a/rb^3/cosB*2*y3b-a/rb^3*
    y2*cotB+3/2*a*y2*W8*cotB/rb^5*2*y3b+y2/rb/W6*(N1*cotB-2*
    nu*y1/W6-a*y1/rb*(1/rb+1/W6))-1/2*y2*W8/rb^3/W6*(N1*
    cotB-2*nu*y1/W6-a*y1/rb*(1/rb+1/W6))*2*y3b-y2*W8/rb/W6^
    2*(N1*cotB-2*nu*y1/W6-a*y1/rb*(1/rb+1/W6))*(y3b/rb+1)+y2*
    W8/rb/W6*(2*nu*y1/W6^2*(y3b/rb+1)+1/2*a*y1/rb^3*(1/rb+
    1/W6)*2*y3b-a*y1/rb*(-1/2/rb^3*2*y3b-1/W6^2*(y3b/rb+
    1)))+y2*cotB/rb/W7*((-2+2*nu)*cosB+W1/W7*W9+a*y3b/rb2/cosB)-
    1/2*y2*W8*cotB/rb^3/W7*((-2+2*nu)*cosB+W1/W7*W9+a*y3b/
    rb2/cosB)*2*y3b-y2*W8*cotB/rb/W7^2*((-2+2*nu)*cosB+W1/W7*
    W9+a*y3b/rb2/cosB)*W3+y2*W8*cotB/rb/W7*((cosB*y3b/rb+1)/
    W7*W9-W1/W7^2*W9*W3-1/2*W1/W7*a/rb^3/cosB*2*y3b+a/rb2/
    cosB-a*y3b/rb2^2/cosB*2*y3b))/pi/(1-nu))+
    b3/2*(1/4*(N1*(-sinB*W3/W7+y1/W6^2*(1+a/rb)*(y3b/rb+1)+
    1/2*y1/W6*a/rb^3*2*y3b+sinB/W7*W2-z1b/W7^2*W2*W3-1/2*
    z1b/W7*a/rb^3*2*y3b)+y1/rb*(a/rb2+1/W6)-1/2*y1*W8/rb^
    3*(a/rb2+1/W6)*2*y3b+y1*W8/rb*(-a/rb2^2*2*y3b-1/W6^2*
    (y3b/rb+1))-1/W7*(sinB*(cosB-a/rb)+z1b/rb*(1+a*y3b/rb2)-1/
    rb/W7*(y2^2*cosB*sinB-a*z1b/rb*W1))+W8/W7^2*(sinB*(cosB-
    a/rb)+z1b/rb*(1+a*y3b/rb2)-1/rb/W7*(y2^2*cosB*sinB-a*z1b/
    rb*W1))*W3-W8/W7*(1/2*sinB*a/rb^3*2*y3b+sinB/rb*(1+a*y3b/
    rb2)-1/2*z1b/rb^3*(1+a*y3b/rb2)*2*y3b+z1b/rb*(a/rb2-a*
    y3b/rb2^2*2*y3b)+1/2/rb^3/W7*(y2^2*cosB*sinB-a*z1b/rb*
    W1)*2*y3b+1/rb/W7^2*(y2^2*cosB*sinB-a*z1b/rb*W1)*W3-1/
    rb/W7*(-a*sinB/rb*W1+1/2*a*z1b/rb^3*W1*2*y3b-a*z1b/rb*
    (cosB*y3b/rb+1))))/pi/(1-nu))+
    b1/2*(1/4*((2-2*nu)*(N1*rFib_ry2*cotB+1/W6*W5-y2^2/W6^2*W5/
    rb-y2^2/W6*a/rb^3-cosB/W7*W2+y2^2*cosB/W7^2*W2/rb+y2^2*
    cosB/W7*a/rb^3)+W8/rb*(2*nu/W6+a/rb2)-y2^2*W8/rb^3*(2*
    nu/W6+a/rb2)+y2*W8/rb*(-2*nu/W6^2/rb*y2-2*a/rb2^2*y2)+
    W8*cosB/rb/W7*(1-2*nu-W1/W7*W2-a*y3b/rb2)-y2^2*W8*cosB/
    rb^3/W7*(1-2*nu-W1/W7*W2-a*y3b/rb2)-y2^2*W8*cosB/rb2/W7^
    2*(1-2*nu-W1/W7*W2-a*y3b/rb2)+y2*W8*cosB/rb/W7*(-1/rb*
    cosB*y2/W7*W2+W1/W7^2*W2/rb*y2+W1/W7*a/rb^3*y2+2*a*
    y3b/rb2^2*y2))/pi/(1-nu))+
    b2/2*(1/4*((-2+2*nu)*N1*cotB*(1/rb*y2/W6-cosB/rb*y2/W7)+(2-
    2*nu)*y1/W6^2*W5/rb*y2+(2-2*nu)*y1/W6*a/rb^3*y2-(2-2*
    nu)*z1b/W7^2*W2/rb*y2-(2-2*nu)*z1b/W7*a/rb^3*y2-W8/rb^
    3*(N1*cotB-2*nu*y1/W6-a*y1/rb2)*y2+W8/rb*(2*nu*y1/W6^2/
    rb*y2+2*a*y1/rb2^2*y2)+W8/W7^2*(cosB*sinB+W1*cotB/rb*((2-
    2*nu)*cosB-W1/W7)+a/rb*(sinB-y3b*z1b/rb2-z1b*W1/rb/W7))/
    rb*y2-W8/W7*(1/rb2*cosB*y2*cotB*((2-2*nu)*cosB-W1/W7)-W1*
    cotB/rb^3*((2-2*nu)*cosB-W1/W7)*y2+W1*cotB/rb*(-cosB/rb*
    y2/W7+W1/W7^2/rb*y2)-a/rb^3*(sinB-y3b*z1b/rb2-z1b*W1/
    rb/W7)*y2+a/rb*(2*y3b*z1b/rb2^2*y2-z1b/rb2*cosB*y2/W7+
    z1b*W1/rb^3/W7*y2+z1b*W1/rb2/W7^2*y2)))/pi/(1-nu))+
    b3/2*(1/4*((2-2*nu)*rFib_ry2+(2-2*nu)*sinB/W7*W2-(2-2*nu)*y2^2*
    sinB/W7^2*W2/rb-(2-2*nu)*y2^2*sinB/W7*a/rb^3+W8*sinB/rb/
    W7*(1+W1/W7*W2+a*y3b/rb2)-y2^2*W8*sinB/rb^3/W7*(1+W1/
    W7*W2+a*y3b/rb2)-y2^2*W8*sinB/rb2/W7^2*(1+W1/W7*W2+a*
    y3b/rb2)+y2*W8*sinB/rb/W7*(1/rb*cosB*y2/W7*W2-W1/W7^2*
    W2/rb*y2-W1/W7*a/rb^3*y2-2*a*y3b/rb2^2*y2))/pi/(1-nu));
	

return(v11[1],v22[1],v33[1],v12[1],v13[1],v23[1])
end