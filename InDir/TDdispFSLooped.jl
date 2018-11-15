function TDdispFSLooped(X,Y,Z,P1,P2,P3,by,bz,bx,nu)
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
eY = [0;1;0];
eZ = [0;0;1];


P1_1=P1[1];P1_2=P1[2];P1_3=P1[3];
P2_1=P2[1];P2_2=P2[2];P2_3=P2[3];
P3_1=P3[1];P3_2=P3[2];P3_3=P3[3];



# Calculate unit strike, dip and normal to TD vectors: For a horizontal TD 
# as an exception, if the normal vector points upward, the strike and dip 
# vectors point Northward and Westward, whereas if the normal vector points
# downward, the strike and dip vectors point Southward and Westward, 
# respectively.
Vnorm = cross(P2-P1,P3-P1);
Vnorm = Vnorm/norm(Vnorm);
Vstrike = cross(eZ,Vnorm);
if norm(Vstrike)==0
	Vstrike = eY*Vnorm[3];
end
Vstrike = Vstrike/norm(Vstrike);
Vdip = cross(Vnorm,Vstrike);

# Transform coordinates from EFCS into TDCS
#At = [Vnorm';Vstrike';Vdip'];
Vx=[Vnorm[1],Vstrike[1],Vdip[1]];
Vy=[Vnorm[2],Vstrike[2],Vdip[2]];
Vz=[Vnorm[3],Vstrike[3],Vdip[3]];
##GARBAGECOLLECTOR STARTS HERE!
(x,y,z)=RotateObject3DNewCoords(X,Y,Z,P2_1,P2_2,P2_3,Vx,Vy,Vz)
(X1,X2,X3)=RotateObject3DNewCoords(P1_1,P1_2,P1_3,P2_1,P2_2,P2_3,Vx,Vy,Vz)
p1=[X1;X2;X3];
(X1,X2,X3)=RotateObject3DNewCoords(P3_1,P3_2,P3_3,P2_1,P2_2,P2_3,Vx,Vy,Vz)
p3=[X1;X2;X3]
p2=zeros(size(p1));
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
A=acos(round(abs(A[1]) ; digits=15));
B=-e12'*e23;
B=acos(round(abs(B[1]) ; digits=15));
C=e23'*e13;
C=acos(round(abs(C[1]) ; digits=15));



#Index out of loop and preallocate
p1yz=p1[2:3];
p2yz=p2[2:3];
p3yz=p3[2:3];
#uV = Array{Float64}(undef, 1,length(X)); 
#vV = Array{Float64}(undef, 1,length(X)); 
#wV = Array{Float64}(undef, 1,length(X)); 

	
me13=-e13;
	
# Determine the best arteact-free configuration for each calculation point
Trimode = trimodefinder(y,z,x,p1yz,p2yz,p3yz);
casepLog=Array{Bool}(undef, length(X),1);
casenLog=Array{Bool}(undef, length(X),1);
casezLog=Array{Bool}(undef, length(X),1);
for i=1:length(Trimode)
	if Trimode[i]==1
		casepLog[i] = true; 
	end
	if Trimode[i]==0;
		casezLog[i] = true;
	end
end


# Calculate first angular dislocation contribution
(u1T,v1T,w1T) = TDSetupD(x,y,z,A,bx,by,bz,nu,p1,me13,casepLog,casenLog,casezLog);
# Calculate second angular dislocation contribution
(u2T,v2T,w2T) = TDSetupD(x,y,z,B,bx,by,bz,nu,p2,e12,casepLog,casenLog,casezLog);
# Calculate third angular dislocation contribution
(u3T,v3T,w3T) = TDSetupD(x,y,z,C,bx,by,bz,nu,p3,e23,casepLog,casenLog,casezLog);		

#=

uelp = AllocSpace(X)
unlp = uelp;
uvlp = uelp;
for i=1:length(X) #For every point in space	

	x=xV[i];
	y=yV[i];
	z=zV[i];
	u=uV[i];
	v=vV[i];
	w=wV[i];
	
	# Calculate the Burgers' function contribution corresponding to the TD
	#Remove indexing

	a = (-x,p1_2-y,p1_3-z);
	b = (-x,-y,-z);
	c = (-x,p3_2-y,p3_3-z);

	na = sqrt(sum(a.^2));
	nb = sqrt(sum(b.^2));
	nc = sqrt(sum(c.^2));
		
	Fi = -2*atan((a[1]*(b[2]*c[3]-b[3]*c[2])-
		a[2]*(b[1]*c[3]-b[3]*c[1])+
		a[3]*(b[1]*c[2]-b[2]*c[1])),
		(na*nb*nc+sum(a.*b)*nc+sum(a.*c)*nb+sum(b.*c)*na))/4/pi;
		
	# Calculate the complete displacement vector components in TDCS
	u = bx*Fi+u;
	v = by*Fi+v;
	w = bz*Fi+w;
	
	# Transform the complete displacement vector components from TDCS into EFCS
	(ue,un,uv) = CoordTrans(u,v,w,At');
	
	uelp[i]=ue;
	unlp[i]=un;
	uvlp[i]=uv;
	
end	

=#
println("X's!")
return(X,X,X)

end


function AllocSpace(X)
Vect = Array{Float64}(undef, length(X),1);
return Vect
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

trimode=Array{Int64}(undef, length(x),1);
for i=1:length(x)

	a = ((p2_2-p3_2)*(x[i]-p3_1)+(p3_1-p2_1)*(y[i]-p3_2))/
		((p2_2-p3_2)*(p1_1-p3_1)+(p3_1-p2_1)*(p1_2-p3_2));
	b = ((p3_2-p1_2)*(x[i]-p3_1)+(p1_1-p3_1)*(y[i]-p3_2))/
		((p2_2-p3_2)*(p1_1-p3_1)+(p3_1-p2_1)*(p1_2-p3_2));
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
	else
		trimode[i]=0;	
	end
		
	if trimode[i]==0 && z!=0
		trimode[i] = 1;
	end

end

return(trimode)
end


function TDSetupD(x,y,z,alpha,bx,by,bz,nu,TriVertex,SideVec,casepLog,casenLog,casezLog)
# TDSetupD transforms coordinates of the calculation points as well as 
# slip vector components from ADCS into TDCS. It then calculates the 
# displacements in ADCS and transforms them into TDCS.



#	if casepLog
#		# Calculate first angular dislocation contribution
#		(u1T,v1T,w1T) = TDSetupD(x,y,z,A,bx,by,bz,nu,p1,me13);
#		# Calculate second angular dislocation contribution
#		(u2T,v2T,w2T) = TDSetupD(x,y,z,B,bx,by,bz,nu,p2,e12);
#		# Calculate third angular dislocation contribution
#		(u3T,v3T,w3T) = TDSetupD(x,y,z,C,bx,by,bz,nu,p3,e23);		
#	elseif casenLog
#		# Calculate first angular dislocation contribution
#		(u1T,v1T,w1T) = TDSetupD(x,y,z,A,bx,by,bz,nu,p1,e13);
#		# Calculate second angular dislocation contribution
#		(u2T,v2T,w2T) = TDSetupD(x,y,z,B,bx,by,bz,nu,p2,me12);
#		# Calculate third angular dislocation contribution
#		(u3T,v3T,w3T) = TDSetupD(x,y,z,C,bx,by,bz,nu,p3,me23);	
#	end


Ct=SideVec[3];
St=SideVec[2];
P1=TriVertex[2];
P2=TriVertex[3];
(y1,z1)  =RotateObject2D(y,z,P1,P2,Ct,St)
(by1,bz1)=RotateObject2D(by,bz,0,0,Ct,St)


if any(casenLog)
	(y1n,z1n)  =RotateObject2D(y,z,P1,P2,Ct,-St) # Unsure this is correct to spin back
	(by1n,bz1n)=RotateObject2D(by,bz,0,0,Ct,-St)
	for i=1:length(x)
		if casenLog[i]==true
			y1[i]=y1n[i];
			z1[i]=z1n[i];
		end
	end
end


### Transformation matrix
#A = [SideVec[3] -SideVec[2];SideVec[2] SideVec[3]];
## Transform coordinates of the calculation points from TDCS into ADCS
#r1 = A*[y-TriVertex[2];z-TriVertex[3]];
#y1 = r1[1];
#z1 = r1[2];
#
## Transform the in-plane slip vector components from TDCS into ADCS
#r2 = A*[by;bz];
#by1 = r2[1];
#bz1 = r2[2];



#Init arrays
Ang=-pi+alpha;
cosA = cos(Ang);
sinA = sin(Ang);
uV  = Array{Float64}(undef, length(x),1);
v0V = Array{Float64}(undef, length(x),1);
w0V = Array{Float64}(undef, length(x),1);

# Calculate displacements associated with an angular dislocation in ADCS
for i=1:length(x)
	if casenLog[i]==true
		(u,v0,w0) = AngDisDisp(x[i],y1[i],z1[i],cosA,sinA,bx,by1n[1],bz1n[1],nu);
	else
		(u,v0,w0) = AngDisDisp(x[i],y1[i],z1[i],cosA,sinA,bx,by1[1],bz1[1],nu);
	end
	uV[i]=u;
	v0V[i]=v0;
	w0V[i]=w0;
end

#=

# Transform displacements from ADCS into TDCS
r3 = A'*[v0;w0];
v = r3[1];
w = r3[2];

# Calculate the "incomplete" displacement vector components in TDCS
if casezLog
	uV[i] = NaN;
	vV[i] = NaN;
	wV[i] = NaN;
else
	uV[i] = u1T+u2T+u3T;
	vV[i] = v1T+v2T+v3T;
	wV[i] = w1T+w2T+w3T;		
end

return(u,v,w)
=#
#println("zeros!")
return(0,0,0)
end


function AngDisDisp(x,y,z,cosA,sinA,bx,by,bz,nu)
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

ux = bx/8/pi/(1-nu)*(x*y/r/(r-z)-x*eta/r/(r-zeta));
vx = bx/8/pi/(1-nu)*(eta*sinA/(r-zeta)-y*eta/r/(r-zeta)+
    y.^2/r/(r-z)+(1-2*nu)*(cosA*log(r-zeta)-log(r-z)));
wx = bx/8/pi/(1-nu)*(eta*cosA/(r-zeta)-y/r-eta*z/r/(r-zeta)-
    (1-2*nu)*sinA*log(r-zeta));
uy = by/8/pi/(1-nu)*(x^2*cosA/r/(r-zeta)-x^2/r/(r-z)-
    (1-2*nu)*(cosA*log(r-zeta)-log(r-z)));
vy = by*x/8/pi/(1-nu)*(y*cosA/r/(r-zeta)-
    sinA*cosA/(r-zeta)-y/r/(r-z));
wy = by*x/8/pi/(1-nu)*(z*cosA/r/(r-zeta)-
    cosA^2/(r-zeta)+1/r);
	
uz = bz*sinA/8/pi/(1-nu)*((1-2*nu)*log(r-zeta)-x.^2/r/(r-zeta));
vz = bz*x*sinA/8/pi/(1-nu)*(sinA/(r-zeta)-y/r/(r-zeta));
wz = bz*x*sinA/8/pi/(1-nu)*(cosA/(r-zeta)-z/r/(r-zeta));

u = ux[1]+uy[1]+uz[1];
v = vx[1]+vy[1]+vz[1];
w = wx[1]+wy[1]+wz[1];

return(u,v,w)
end
