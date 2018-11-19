function TDdispFS(X,Y,Z,P1,P2,P3,by,bz,bx,nu)
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
# ue, un and uv:
# Calculated displacement vector components in EFCS. ue, un and uv have
# the same unit as Ss, Ds and Ts in the inputs.
# 
# 
# Example: Calculate and plot the first component of displacement vector 
# on a regular grid.
# 
# (X,Y,Z) = meshgrid(-3:.02:3,-3:.02:3,2);
# (ue,un,uv) = TDdispFS(X,Y,Z,(-1 0 0),(1 -1 -1),(0 1.5 .5),-1,2,3,.25);
# h = surf(X,Y,reshape(ue,size(X)),'edgecolor','none');
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

(x,y,z)=RotateObject3DNewCoords(X,Y,Z,P2_1,P2_2,P2_3,Vx,Vy,Vz)
(X1,X2,X3)=RotateObject3DNewCoords(P1_1,P1_2,P1_3,P2_1,P2_2,P2_3,Vx,Vy,Vz)
p1=[X1;X2;X3];
(X1,X2,X3)=RotateObject3DNewCoords(P3_1,P3_2,P3_3,P2_1,P2_2,P2_3,Vx,Vy,Vz)
p3=[X1;X2;X3]
p1_2=p1[2];p1_3=p1[3];
p3_2=p3[2];p3_3=p3[3];


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
	
# Determine the best arteact-free configuration for each calculation point
Trimode = trimodefinder(y,z,x,p1[2:3],p2[2:3],p3[2:3]);
casepLog=falses(length(X),1);
casenLog=falses(length(X),1);
casezLog=falses(length(X),1);
for i=1:length(Trimode)
	if Trimode[i]==1
		casepLog[i] = 1; 
	end
	if Trimode[i]==-1
		casenLog[i] = 1; 
	end	
	if Trimode[i]==0;
		casezLog[i] = 1;
	end
end

#Assuming Configuration I
me13=-e13;
# Calculate first angular dislocation contribution
(u1T,v1T,w1T) = TDSetupD(x,y,z,A,bx,by,bz,nu,p1,me13,casenLog);
# Calculate second angular dislocation contribution
(u2T,v2T,w2T) = TDSetupD(x,y,z,B,bx,by,bz,nu,p2,e12 ,casenLog);
# Calculate third angular dislocation contribution
(u3T,v3T,w3T) = TDSetupD(x,y,z,C,bx,by,bz,nu,p3,e23 ,casenLog);		

#Do some allocation before loop
#a1=-x;	a2=p1_2.-y;	a3=p1_3.-z;
#b1=-x;	b2=-y;		b3=-z;
#c1=-x;	c2=p3_2.-y;	c3=p3_3.-z;
uV = Array{Float64}(undef, length(X),1); 
vV = Array{Float64}(undef, length(X),1); 
wV = Array{Float64}(undef, length(X),1); 

# Calculate the "incomplete" displacement vector components in TDCS
for i=1:length(x)
	if casezLog[i] == 1; 
		uV[i] = NaN;
		vV[i] = NaN;
		wV[i] = NaN;
	else
		uV[i] = u1T[i]+u2T[i]+u3T[i];
		vV[i] = v1T[i]+v2T[i]+v3T[i];
		wV[i] = w1T[i]+w2T[i]+w3T[i];
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
	uV[i] = bx*Fi+uV[i];
	vV[i] = by*Fi+vV[i];
	wV[i] = bz*Fi+wV[i];
	
end	

# Transform the complete displacement vector components from TDCS into EFCS
(ue,un,uv)=RotateObject3DNewCoords(uV,vV,wV,0,0,0,Vnorm,Vstrike,Vdip)

return(ue,un,uv)

end






