function TDstrainFS(X,Y,Z,P1,P2,P3,by,bz,bx,mu,lambda)
# TDstressFS 
# Calculates stresses and strains associated with a triangular dislocation 
# in an elastic full-space.
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
# mu and lambda:
# Lame constants.
#
# OUTPUTS
# Stress:
# Calculated stress tensor components in EFCS. The six columns of Stress 
# are Sxx, Syy, Szz, Sxy, Sxz and Syz, respectively. The stress components 
# have the same unit as Lame constants.
#
# Strain:
# Calculated strain tensor components in EFCS. The six columns of Strain 
# are Exx, Eyy, Ezz, Exy, Exz and Eyz, respectively. The strain components 
# are dimensionless.
# 
# 
# Example: Calculate and plot the first component of stress tensor on a  
# regular grid.
# 
# [X,Y,Z] = meshgrid(-3:.02:3,-3:.02:3,2);
# [Stress,Strain] = TDstressFS(X,Y,Z,[-1 0 0],[1 -1 -1],[0 1.5 .5],...
# -1,2,3,.33e11,.33e11);
# h = surf(X,Y,reshape(Stress(:,1),size(X)),'edgecolor','none');
# view(2)
# axis equal
# axis tight
# set(gcf,'renderer','painters')

# Reference journal article: 
# Nikkhoo M. and Walter T.R., 2015. Triangular dislocation: An analytical, 
# artefact-free solution. 
# Submitted to Geophysical Journal International 

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

# I appreciate any comments or bug reports.

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
# Email: 
# mehdi.nikkhoo@gfz-potsdam.de 
# mehdi.nikkhoo@gmail.com

nu =lambda/(2*(mu+lambda));    #Poisson's ratio, Equation 8.28 Pollard

#bx = Ts; # Tensile-slip
#by = Ss; # Strike-slip
#bz = Ds; # Dip-slip

#Init some vars
p1 = zeros(3,1);
p2 = zeros(3,1);
p3 = zeros(3,1);
P1_1=P1[1];P1_2=P1[2];P1_3=P1[3];
P2_1=P2[1];P2_2=P2[2];P2_3=P2[3];
P3_1=P3[1];P3_2=P3[2];P3_3=P3[3];

(Vnorm,Vstrike,Vdip)=CalculateLocalTriCoords(P1,P2,P3)

# Transform coordinates from EFCS into TDCS
#At = [Vnorm';Vstrike';Vdip'];
Vx=[Vnorm[1],Vstrike[1],Vdip[1]];
Vy=[Vnorm[2],Vstrike[2],Vdip[2]];
Vz=[Vnorm[3],Vstrike[3],Vdip[3]];


# Transform coordinates from EFCS into TDCS
#At = [Vnorm';Vstrike';Vdip'];
Vx=[Vnorm[1],Vstrike[1],Vdip[1]];
Vy=[Vnorm[2],Vstrike[2],Vdip[2]];
Vz=[Vnorm[3],Vstrike[3],Vdip[3]];

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
@info "Some variables" X[90],Y[90]
poop
(casepLog,casenLog,casezLog) = trimodefinder(y,z,x,p1[2:3],p2[2:3],p3[2:3]);

xp=x[casepLog];
yp=y[casepLog];
zp=z[casepLog];
xn=x[casenLog];
yn=y[casenLog];
zn=z[casenLog];




# Calculate first angular dislocation contribution POS
(Exx1Tp,Eyy1Tp,Ezz1Tp,Exy1Tp,Exz1Tp,Eyz1Tp) = TDSetupS(xp,yp,zp,A,bx,by,bz,nu,p1,-e13);
# Calculate second angular dislocation contribution
(Exx2Tp,Eyy2Tp,Ezz2Tp,Exy2Tp,Exz2Tp,Eyz2Tp) = TDSetupS(xp,yp,zp,B,bx,by,bz,nu,p2,e12);
# Calculate third angular dislocation contribution
(Exx3Tp,Eyy3Tp,Ezz3Tp,Exy3Tp,Exz3Tp,Eyz3Tp) = TDSetupS(xp,yp,zp,C,bx,by,bz,nu,p3,e23);

# Calculate first angular dislocation contribution NEG
(Exx1Tn,Eyy1Tn,Ezz1Tn,Exy1Tn,Exz1Tn,Eyz1Tn) = TDSetupS(xn,yn,zn,A,bx,by,bz,nu,p1,e13);
# Calculate second angular dislocation contribution
(Exx2Tn,Eyy2Tn,Ezz2Tn,Exy2Tn,Exz2Tn,Eyz2Tn) = TDSetupS(xn,yn,zn,B,bx,by,bz,nu,p2,-e12);
# Calculate third angular dislocation contribution
(Exx3Tn,Eyy3Tn,Ezz3Tn,Exy3Tn,Exz3Tn,Eyz3Tn) = TDSetupS(xn,yn,zn,C,bx,by,bz,nu,p3,-e23);	

#Do some allocation before loop
exx = Array{Float64}(undef, length(X),1); 
eyy = Array{Float64}(undef, length(X),1); 
ezz = Array{Float64}(undef, length(X),1); 
exy = Array{Float64}(undef, length(X),1); 
exz = Array{Float64}(undef, length(X),1); 
eyz = Array{Float64}(undef, length(X),1); 
		
exx[casenLog]=Exx1Tn.+Exx2Tn.+Exx3Tn;
eyy[casenLog]=Eyy1Tn.+Eyy2Tn.+Eyy3Tn;
ezz[casenLog]=Ezz1Tn.+Ezz2Tn.+Ezz3Tn;
exy[casenLog]=Exy1Tn.+Exy2Tn.+Exy3Tn;
exz[casenLog]=Exz1Tn.+Exz2Tn.+Exz3Tn;
eyz[casenLog]=Eyz1Tn.+Eyz2Tn.+Eyz3Tn;

exx[casepLog]=Exx1Tp.+Exx2Tp.+Exx3Tp;
eyy[casepLog]=Eyy1Tp.+Eyy2Tp.+Eyy3Tp;
ezz[casepLog]=Ezz1Tp.+Ezz2Tp.+Ezz3Tp;
exy[casepLog]=Exy1Tp.+Exy2Tp.+Exy3Tp;
exz[casepLog]=Exz1Tp.+Exz2Tp.+Exz3Tp;
eyz[casepLog]=Eyz1Tp.+Eyz2Tp.+Eyz3Tp;	
		
# Calculate the strain tensor components in TDCS
for i=1:length(x)
	if casezLog[i] == 1; 
		exx[i] = NaN;
		eyy[i] = NaN;
		ezz[i] = NaN;
		exy[i] = NaN;
		exz[i] = NaN;
		eyz[i] = NaN;
	end
end

# Transform the strain tensor components from TDCS into EFCS
(Exx,Eyy,Ezz,Exy,Exz,Eyz) = TensorTransformation3D(exx,eyy,ezz,exy,exz,eyz,[Vnorm Vstrike Vdip]);


return(Exx,Eyy,Ezz,Exy,Exz,Eyz)

end