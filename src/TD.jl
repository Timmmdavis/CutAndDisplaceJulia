function TD(X::Array{Float64,2},Y::Array{Float64,2},Z::Array{Float64,2},
		    P1List::Array{Float64,2},P2List::Array{Float64,2},P3List::Array{Float64,2},
			Dss::Array,Dds::Array,Dn::Array,nu::Float64,mu::Float64,
			DispFlag::Int64,StrainFlag::Int64,HSflag::Int64)
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

###println("This initial allocation should follow through the functions....")
#Do some allocation before loop
if DispFlag==1
	
	#First (and hopefully last allocation)
	UxDn = zeros(length(X),SzCmp); 
	UyDn = zeros(length(X),SzCmp); 
	UzDn = zeros(length(X),SzCmp); 
	UxDss = zeros(length(X),SzCmp); 
	UyDss = zeros(length(X),SzCmp); 
	UzDss = zeros(length(X),SzCmp); 
	UxDds = zeros(length(X),SzCmp); 
	UyDds = zeros(length(X),SzCmp); 
	UzDds = zeros(length(X),SzCmp); 
	
	zers= zeros(length(X),1); 
	
end
if StrainFlag==1

	#First (and hopefully last allocation)
	ExxDn = zeros(length(X),SzCmp); 
	EyyDn = zeros(length(X),SzCmp); 
	EzzDn = zeros(length(X),SzCmp);
	ExyDn = zeros(length(X),SzCmp); 
	ExzDn = zeros(length(X),SzCmp); 
	EyzDn = zeros(length(X),SzCmp);
	ExxDss = zeros(length(X),SzCmp); 
	EyyDss = zeros(length(X),SzCmp); 
	EzzDss = zeros(length(X),SzCmp);
	ExyDss = zeros(length(X),SzCmp); 
	ExzDss = zeros(length(X),SzCmp); 
	EyzDss = zeros(length(X),SzCmp);
	ExxDds = zeros(length(X),SzCmp); 
	EyyDds = zeros(length(X),SzCmp); 
	EzzDds = zeros(length(X),SzCmp);
	ExyDds = zeros(length(X),SzCmp); 
	ExzDds = zeros(length(X),SzCmp); 
	EyzDds = zeros(length(X),SzCmp);
end

for i=1:SzCmp #For every element (multithreaded)  Threads.@threads 

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
		P1i[3] = copy(-P1[3]);
		P2i[3] = copy(-P2[3]);
		P3i[3] = copy(-P3[3]);
	end

	if DispFlag==1

		# Calculate main dislocation contribution to displacements
		(UxDn[:,i], UyDn[:,i], UzDn[:,i],
		 UxDss[:,i],UyDss[:,i],UzDss[:,i],
		 UxDds[:,i],UyDds[:,i],UzDds[:,i]) = 
		 TDdispFS(X,Y,Z,P1,P2,P3,Dss[i],Dds[i],Dn[i],nu,0,
		 UxDn[:,i], UyDn[:,i], UzDn[:,i],
		 UxDss[:,i],UyDss[:,i],UzDss[:,i],
		 UxDds[:,i],UyDds[:,i],UzDds[:,i]);
		
		if HSflag==1
			
			
			# Calculate harmonic fUyction contribution to displacements
			(UxDn[:,i], UyDn[:,i], UzDn[:,i],
			 UxDss[:,i],UyDss[:,i],UzDss[:,i],
		     UxDds[:,i],UyDds[:,i],UzDds[:,i]) = 
			 TDdisp_HarFunc(X,Y,Z,P1,P2,P3,Dss[i],Dds[i],Dn[i],nu,
			 UxDn[:,i], UyDn[:,i], UzDn[:,i],
			 UxDss[:,i],UyDss[:,i],UzDss[:,i],
			 UxDds[:,i],UyDds[:,i],UzDds[:,i]);
			
			# Calculate image dislocation contribution to displacements
			(UxDn[:,i], UyDn[:,i], UzDn[:,i],
			 UxDss[:,i],UyDss[:,i],UzDss[:,i],
		     UxDds[:,i],UyDds[:,i],UzDds[:,i]) =
			 TDdispFS(X,Y,Z,P1i,P2i,P3i,Dss[i],Dds[i],Dn[i],nu,1,
			 UxDn[:,i], UyDn[:,i], UzDn[:,i],
			 UxDss[:,i],UyDss[:,i],UzDss[:,i],
			 UxDds[:,i],UyDds[:,i],UzDds[:,i]);
					
		end #HsFlag

	end#DispFlag 

	if StrainFlag==1

		#Elastic con
		lambda=(2*mu*nu)/(1-(2*nu));
		
		(ExxDn[:,i], EyyDn[:,i], EzzDn[:,i], ExyDn[:,i], ExzDn[:,i], EyzDn[:,i],
		 ExxDss[:,i],EyyDss[:,i],EzzDss[:,i],ExyDss[:,i],ExzDss[:,i],EyzDss[:,i],
		 ExxDds[:,i],EyyDds[:,i],EzzDds[:,i],ExyDds[:,i],ExzDds[:,i],EyzDds[:,i]) =
		 TDstrainFS(X,Y,Z,P1,P2,P3,Dss[i],Dds[i],Dn[i],mu,lambda,0,
		 ExxDn[:,i], EyyDn[:,i], EzzDn[:,i], ExyDn[:,i], ExzDn[:,i], EyzDn[:,i],
		 ExxDss[:,i],EyyDss[:,i],EzzDss[:,i],ExyDss[:,i],ExzDss[:,i],EyzDss[:,i],
		 ExxDds[:,i],EyyDds[:,i],EzzDds[:,i],ExyDds[:,i],ExzDds[:,i],EyzDds[:,i]);	
		
		if HSflag==1
		
			(ExxDn[:,i], EyyDn[:,i], EzzDn[:,i], ExyDn[:,i], ExzDn[:,i], EyzDn[:,i],
		     ExxDss[:,i],EyyDss[:,i],EzzDss[:,i],ExyDss[:,i],ExzDss[:,i],EyzDss[:,i],
		     ExxDds[:,i],EyyDds[:,i],EzzDds[:,i],ExyDds[:,i],ExzDds[:,i],EyzDds[:,i]) =
			 TDstrain_HarFunc(X,Y,Z,P1,P2,P3,Dss[i],Dds[i],Dn[i],mu,lambda,nu,
			 ExxDn[:,i], EyyDn[:,i], EzzDn[:,i], ExyDn[:,i], ExzDn[:,i], EyzDn[:,i],
			 ExxDss[:,i],EyyDss[:,i],EzzDss[:,i],ExyDss[:,i],ExzDss[:,i],EyzDss[:,i],
			 ExxDds[:,i],EyyDds[:,i],EzzDds[:,i],ExyDds[:,i],ExzDds[:,i],EyzDds[:,i]);	
		
			# Calculate image dislocation contribution to strains and stresses
			(ExxDn[:,i], EyyDn[:,i], EzzDn[:,i], ExyDn[:,i], ExzDn[:,i], EyzDn[:,i],
			 ExxDss[:,i],EyyDss[:,i],EzzDss[:,i],ExyDss[:,i],ExzDss[:,i],EyzDss[:,i],
			 ExxDds[:,i],EyyDds[:,i],EzzDds[:,i],ExyDds[:,i],ExzDds[:,i],EyzDds[:,i]) =
			 TDstrainFS(X,Y,Z,P1i,P2i,P3i,Dss[i],Dds[i],Dn[i],mu,lambda,1,
			 ExxDn[:,i], EyyDn[:,i], EzzDn[:,i], ExyDn[:,i], ExzDn[:,i], EyzDn[:,i],
			 ExxDss[:,i],EyyDss[:,i],EzzDss[:,i],ExyDss[:,i],ExzDss[:,i],EyzDss[:,i],
			 ExxDds[:,i],EyyDds[:,i],EzzDds[:,i],ExyDds[:,i],ExzDds[:,i],EyzDds[:,i]);
			
		end #HS flag
		
	end #Stress flag
	
end #Over all els
	
	
if DispFlag==0
	return(ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
		   ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
		   ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds)
elseif StrainFlag==0
	return(UxDn,UyDn,UzDn,
		   UxDss,UyDss,UzDss,
		   UxDds,UyDds,UzDds)
else

	return(ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
		   ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
		   ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,
		   UxDn,UyDn,UzDn,
		   UxDss,UyDss,UzDss,
		   UxDds,UyDds,UzDds)
	
end
end


function TDdispFS(X,Y,Z,P1,P2,P3,Dss,Dds,Dn,nu,ImageFlag,
				UxDn, UyDn, UzDn,
				UxDss,UyDss,UzDss,
				UxDds,UyDds,UzDds)
#Calculate displacement components in a full space
(Vnorm,Vstrike,Vdip)=CalculateLocalTriCoords(P1,P2,P3)
(p1,p2,p3,x,y,z,e12,e13,e23,A,B,C,Pos,Neg,casezLog)=GlobalToTDECoords(P1,P2,P3,X,Y,Z,Vnorm,Vstrike,Vdip)

#Only do on 2nd time we go into func, there is no point rotating 0's centred at 0...
if ImageFlag==1; #This means we are computing the iamge dislocation
	#Inverse rot mat
	VxR=[Vnorm[1],Vstrike[1],Vdip[1]];
	VyR=[Vnorm[2],Vstrike[2],Vdip[2]];
	VzR=[Vnorm[3],Vstrike[3],Vdip[3]];
	# Transform the complete displacement vector EFCS to TDCS 
	(UxDn, UyDn, UzDn) =RotateObject3DNewCoords!(UxDn, UyDn, UzDn ,0,0,0,VxR,VyR,VzR)
	(UxDss,UyDss,UzDss)=RotateObject3DNewCoords!(UxDss,UyDss,UzDss,0,0,0,VxR,VyR,VzR)
	(UxDds,UyDds,UzDds)=RotateObject3DNewCoords!(UxDds,UyDds,UzDds,0,0,0,VxR,VyR,VzR)
	
	#Flip these so we add component
	if P1[3]==0 && P2[3]==0 && P3[3]==0
		UzDn  = -UzDn;
		UzDss = -UzDss;
		UzDds = -UzDds;
	end
end


#Turn to cart index
Pos=findall(Pos);
Neg=findall(Neg)

# Calculate first angular dislocation contribution POS
(UxDn[Pos],UxDss[Pos],UxDds[Pos],
 UyDn[Pos],UyDss[Pos],UyDds[Pos],
 UzDn[Pos],UzDss[Pos],UzDds[Pos]) = 
 TDSetupD(x[Pos],y[Pos],z[Pos],A,Dn,Dss,Dds,nu,p1,-e13,
 UxDn[Pos],UxDss[Pos],UxDds[Pos],
 UyDn[Pos],UyDss[Pos],UyDds[Pos],
 UzDn[Pos],UzDss[Pos],UzDds[Pos]);
# Calculate second angular dislocation contribution
(UxDn[Pos],UxDss[Pos],UxDds[Pos],
 UyDn[Pos],UyDss[Pos],UyDds[Pos],
 UzDn[Pos],UzDss[Pos],UzDds[Pos]) =
 TDSetupD(x[Pos],y[Pos],z[Pos],B,Dn,Dss,Dds,nu,p2,e12,
 UxDn[Pos],UxDss[Pos],UxDds[Pos],
 UyDn[Pos],UyDss[Pos],UyDds[Pos],
 UzDn[Pos],UzDss[Pos],UzDds[Pos]); 
# Calculate third angular dislocation contribution
(UxDn[Pos],UxDss[Pos],UxDds[Pos],
 UyDn[Pos],UyDss[Pos],UyDds[Pos],
 UzDn[Pos],UzDss[Pos],UzDds[Pos]) =
 TDSetupD(x[Pos],y[Pos],z[Pos],C,Dn,Dss,Dds,nu,p3,e23,
 UxDn[Pos],UxDss[Pos],UxDds[Pos],
 UyDn[Pos],UyDss[Pos],UyDds[Pos],
 UzDn[Pos],UzDss[Pos],UzDds[Pos]);

 
# Calculate first angular dislocation contribution NEG
(UxDn[Neg],UxDss[Neg],UxDds[Neg],
 UyDn[Neg],UyDss[Neg],UyDds[Neg],
 UzDn[Neg],UzDss[Neg],UzDds[Neg]) =
 TDSetupD(x[Neg],y[Neg],z[Neg],A,Dn,Dss,Dds,nu,p1,e13,
 UxDn[Neg],UxDss[Neg],UxDds[Neg],
 UyDn[Neg],UyDss[Neg],UyDds[Neg],
 UzDn[Neg],UzDss[Neg],UzDds[Neg]);
# Calculate second angular dislocation contribution
(UxDn[Neg],UxDss[Neg],UxDds[Neg],
 UyDn[Neg],UyDss[Neg],UyDds[Neg],
 UzDn[Neg],UzDss[Neg],UzDds[Neg]) = 
 TDSetupD(x[Neg],y[Neg],z[Neg],B,Dn,Dss,Dds,nu,p2,-e12,
 UxDn[Neg],UxDss[Neg],UxDds[Neg],
 UyDn[Neg],UyDss[Neg],UyDds[Neg],
 UzDn[Neg],UzDss[Neg],UzDds[Neg]);
# Calculate third angular dislocation contribution
(UxDn[Neg],UxDss[Neg],UxDds[Neg],
 UyDn[Neg],UyDss[Neg],UyDds[Neg],
 UzDn[Neg],UzDss[Neg],UzDds[Neg]) =
 TDSetupD(x[Neg],y[Neg],z[Neg],C,Dn,Dss,Dds,nu,p3,-e23,
 UxDn[Neg],UxDss[Neg],UxDds[Neg],
 UyDn[Neg],UyDss[Neg],UyDds[Neg],
 UzDn[Neg],UzDss[Neg],UzDds[Neg]);	
 

# Calculate the "incomplete" displacement vector components in TDCS
for i=1:length(x)

	if casezLog[i] == 1; 
		UxDn[i] = NaN;
		UyDn[i] = NaN;
		UzDn[i] = NaN;
		UxDss[i] = NaN;
		UyDss[i] = NaN;
		UzDss[i] = NaN;
		UxDds[i] = NaN;
		UyDds[i] = NaN;
		UzDds[i] = NaN;
	end

	a1=-x[i];	a2=p1[2]-y[i];	a3=p1[3]-z[i];
	b1=-x[i];	b2=-y[i];		b3=-z[i];
	c1=-x[i];	c2=p3[2]-y[i];	c3=p3[3]-z[i];
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
	UxDn[i]  = Dn*Fi+UxDn[i];
	UyDss[i] = Dss*Fi+UyDss[i];
	UzDds[i] = Dds*Fi+UzDds[i];
	#Only add parts that matter

	
end	


# Transform the complete displacement vector components from TDCS into EFCS
(UxDn, UyDn, UzDn) =RotateObject3DNewCoords!(UxDn, UyDn, UzDn ,0,0,0,Vnorm,Vstrike,Vdip)
(UxDss,UyDss,UzDss)=RotateObject3DNewCoords!(UxDss,UyDss,UzDss,0,0,0,Vnorm,Vstrike,Vdip)
(UxDds,UyDds,UzDds)=RotateObject3DNewCoords!(UxDds,UyDds,UzDds,0,0,0,Vnorm,Vstrike,Vdip)

#Flip sign if flat tri 
if ImageFlag==1; #This means we are computing the iamge dislocation	
	if P1[3]==0 && P2[3]==0 && P3[3]==0
		UxDn  = -UxDn;
		UyDn  = -UyDn;
		UzDn  = -UzDn;
		UxDss = -UxDss;
		UyDss = -UyDss;
		UzDss = -UzDss;
		UxDds = -UxDds;
		UyDds = -UyDds;		
		UzDds = -UzDds;		
	end
end

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

(x,y,z)=RotateObject3DNewCoords(X,Y,Z,P2[1],P2[2],P2[3],Vx,Vy,Vz)
(p1[1],p1[2],p1[3])=RotateObject3DNewCoords(P1[1],P1[2],P1[3],P2[1],P2[2],P2[3],Vx,Vy,Vz)
(p3[1],p3[2],p3[3])=RotateObject3DNewCoords(P3[1],P3[2],P3[3],P2[1],P2[2],P2[3],Vx,Vy,Vz)#Vx,Vy,Vz

#Get interior angles and vectors along the triangle edges. 
(e12,e13,e23,A,B,C)=CalcTDVectsAndAngles(p1,p2,p3)
	
# Determine the best arteact-free configuration for each calculation point
(casepLog,casenLog,casezLog) = trimodefinder(y,z,x,p1[2:3],p2[2:3],p3[2:3]);

return(p1,p2,p3,x,y,z,e12,e13,e23,A,B,C,casepLog,casenLog,casezLog)
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



function TDSetupD(x,y,z,alpha,Dn,Dss,Dds,nu,TriVertex,SideVec,
 UxDn,UxDss,UxDds,
 UyDn,UyDss,UyDds,
 UzDn,UzDss,UzDds)
# TDSetupD transforms coordinates of the calculation points as well as 
# slip vector components from ADCS into TDCS. It then calculates the 
# displacements in ADCS and transforms them into TDCS.

Ct=SideVec[3];
St=SideVec[2];

# Transform coordinates of the calculation points from TDCS into ADCS
(y,z)  =RotateObject2D!(y,z,TriVertex[2],TriVertex[3],Ct,St)
# Transform the in-plane slip vector components from TDCS into ADCS
(Dss1,Dds0)=RotateObject2D(Dss,0.,0,0,Ct,St)
(Dss0,Dds1)=RotateObject2D(0.,Dds,0,0,Ct,St)


#Init arrays
Ang=-pi+alpha;
cosA = cos(Ang);
sinA = sin(Ang);

#Extra defs out of loop to speed it up
E1=(1-nu); #Elastic cons
E2=(1-2*nu);
cosA2=cosA^2;
sinADE1=sinA/8/pi/(1-nu);

Dn8p=Dn/8/pi;

# Calculate displacements associated with an angular dislocation in ADCS
for i=1:length(x)
	
	(ux,uy,uz,vx,vy,vz,wx,wy,wz) = AngDisDisp(x[i],y[i],z[i],cosA,sinA,E1,E2,cosA2,sinADE1);

	#Ux is in global coords:
	#components due to opening
	UxDn[i]=UxDn[i]+(Dn8p/E1*ux);
	#Comp mixed components, Ux is in global coords too
	UxDss[i]=UxDss[i]+(Dss1/8/pi/E1*uy)+(Dds0*sinADE1*uz)
	UxDds[i]=UxDds[i]+(Dss0/8/pi/E1*uy)+(Dds1*sinADE1*uz)	
	#Comp Uy and Uz in the current coords and rotate back
	uyDn=(Dn8p/E1*vx);
	uzDn=(Dn8p/E1*wx);
	#Comp mixed components 
	uyDss=(Dss1*x[i]/8/pi/E1*vy)+(Dds0*x[i]*sinADE1*vz)	
	uyDds=(Dss0*x[i]/8/pi/E1*vy)+(Dds1*x[i]*sinADE1*vz)		
	uzDss=(Dss1*x[i]/8/pi/E1*wy)+(Dds0*x[i]*sinADE1*wz)
	uzDds=(Dss0*x[i]/8/pi/E1*wy)+(Dds1*x[i]*sinADE1*wz)
	#Rotate these to global coords		
	(uyDn, uzDn)= RotateObject2D!(uyDn, uzDn ,0,0,Ct,-St)	
	(uyDss,uzDss)=RotateObject2D!(uyDss,uzDss,0,0,Ct,-St) #Rotate back
	(uyDds,uzDds)=RotateObject2D!(uyDds,uzDds,0,0,Ct,-St) #Rotate back
	#Add these to our total vector 
	UyDn[i]=UyDn[i]+uyDn;
	UzDn[i]=UzDn[i]+uzDn;
	UyDss[i]=UyDss[i]+uyDss;
	UyDds[i]=UyDds[i]+uyDds;
	UzDss[i]=UzDss[i]+uzDss;
	UzDds[i]=UzDds[i]+uzDds;

	
end

## Transform coordinates of the calculation points from ADCS into TDCS
(y,z)  =RotateObject2D!(y,z,0,0,Ct,-St)
y.=y.-TriVertex[2];
z.=z.-TriVertex[3];

return( UxDn,UxDss,UxDds,
		UyDn,UyDss,UyDds,
		UzDn,UzDss,UzDds)
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

rMzeta=r-zeta;
rMz=r-z;

ux = x*y/r/rMz-x*eta/r/rMzeta; 													#b8p/E1*
vx = eta*sinA/rMzeta-y*eta/r/rMzeta+y.^2/r/rMz+E2*(cosA*log(rMzeta)-log(rMz));	#b8p/E1*
wx = eta*cosA/rMzeta-y/r-eta*z/r/rMzeta-E2*sinA*log(rMzeta);					#b8p/E1*
	
uy = x^2*cosA/r/rMzeta-x^2/r/rMz-E2*(cosA*log(rMzeta)-log(rMz)); 				#by/8/pi/E1*
vy = y*cosA/r/rMzeta-sinA*cosA/rMzeta-y/r/rMz;									#by*x/8/pi/E1*			
wy = z*cosA/r/rMzeta-cosA2/rMzeta+1/r;											#by*x/8/pi/E1*
	
uz = E2*log(rMzeta)-x.^2/r/rMzeta;												#bz*sinADE1*
vz = sinA/rMzeta-y/r/rMzeta;													#bz*x*sinADE1*	
wz = cosA/rMzeta-z/r/rMzeta;													#bz*x*sinADE1*


#Export individual components
return(ux[1],uy[1],uz[1],vx[1],vy[1],vz[1],wx[1],wy[1],wz[1])
end

function TDdisp_HarFunc(X,Y,Z,P1,P2,P3,Dss,Dds,Dn,nu,
				UxDn, UyDn, UzDn,
				UxDss,UyDss,UzDss,
				UxDds,UyDds,UzDds)
# TDdisp_HarFunc calculates the harmonic function contribution to the
# displacements associated with a triangular dislocation in a half-space.
# The function cancels the surface normal tractions induced by the main and
# image dislocations.

# Calculate unit strike, dip and normal to TD vectors: 
(Vnorm,Vstrike,Vdip)=CalculateLocalTriCoords(P1,P2,P3)

# Calculate contribution of angular dislocation pair on each TD side 

# Side P1P2
(UxDn, UyDn, UzDn,
 UxDss,UyDss,UzDss,
 UxDds,UyDds,UzDds) =
 AngSetupDispFSC(X,Y,Z,P1,P2,nu,Vnorm,Vstrike,Vdip,Dn,Dss,Dds,
 UxDn, UyDn, UzDn,
 UxDss,UyDss,UzDss,
 UxDds,UyDds,UzDds); 
 
# Side P2P3
(UxDn, UyDn, UzDn,
 UxDss,UyDss,UzDss,
 UxDds,UyDds,UzDds) = 
 AngSetupDispFSC(X,Y,Z,P2,P3,nu,Vnorm,Vstrike,Vdip,Dn,Dss,Dds,
 UxDn, UyDn, UzDn,
 UxDss,UyDss,UzDss,
 UxDds,UyDds,UzDds); 
 
# Side P3P1
(UxDn, UyDn, UzDn,
 UxDss,UyDss,UzDss,
 UxDds,UyDds,UzDds) =
 AngSetupDispFSC(X,Y,Z,P3,P1,nu,Vnorm,Vstrike,Vdip,Dn,Dss,Dds,
 UxDn, UyDn, UzDn,
 UxDss,UyDss,UzDss,
 UxDds,UyDds,UzDds); 

return(UxDn,UyDn,UzDn,UxDss,UyDss,UzDss,UxDds,UyDds,UzDds);
end


function AngSetupDispFSC(X,Y,Z,PA,PB,nu,Vnorm,Vstrike,Vdip,Dn,Dss,Dds,
 UxDn,UyDn,UzDn,
 UxDss,UyDss,UzDss,
 UxDds,UyDds,UzDds)
# AngSetupFSC calculates the Free Surface Correction to displacements 
# associated with angular dislocation pair on each TD side.

# Calculate TD side vector and the angle of the angular dislocation pair
(SideVec,eZ,beta)=CalcSideVec(PA,PB)

if abs(beta)<eps() || abs(pi-beta)<eps()
	#Simply add nothing to the values coming in
else
    (I,y1A,y2A,y3A,y1B,y2B,y3B,ey1,ey2,ey3)=CalcSlipVectorDiscCoords(SideVec,eZ,X,Y,Z,PA,beta)
    
	#Inverse rot mat
	VxR=[ey1[1],ey2[1],ey3[1]];
	VyR=[ey1[2],ey2[2],ey3[2]];
	VzR=[ey1[3],ey2[3],ey3[3]];

	## Transform slip vector components from TDCS into EFCS
	(Dn1__,Dss0n_,Dds0n_) = RotateObject3DNewCoords(Dn,0.,0.,0,0,0,Vnorm,Vstrike,Vdip);
	(Dn0ss,Dss1__,Dds1ss) = RotateObject3DNewCoords(0.,Dss,0.,0,0,0,Vnorm,Vstrike,Vdip);
	(Dn0ds,Dss0ds,Dds1__) = RotateObject3DNewCoords(0.,0.,Dds,0,0,0,Vnorm,Vstrike,Vdip);
	# Transform slip vector components from EFCS to ADCS
	(Dn1__,Dss0n_,Dds0n_)=RotateObject3DNewCoords(Dn1__,Dss0n_,Dds0n_,0,0,0,ey1,ey2,ey3)
	(Dn0ss,Dss1__,Dds1ss)=RotateObject3DNewCoords(Dn0ss,Dss1__,Dds1ss,0,0,0,ey1,ey2,ey3)
	(Dn0ds,Dss0ds,Dds1__)=RotateObject3DNewCoords(Dn0ds,Dss0ds,Dds1__,0,0,0,ey1,ey2,ey3)

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
	
		(uxA,uyA,uzA,vxA,vyA,vzA,wxA,wyA,wzA) = AngDisDispFSC(y1A[indx[i]],y2A[indx[i]],y3A[indx[i]],cosB,sinB,cotB,cotB2,nu,-PA[3]);
		
		
		#Add mixed components together (different coords)
		uxDn =-((Dn1__/4/pi/(1-nu)*uxA)+(Dss0n_/4/pi/(1-nu)*uyA)+(Dds0n_/4/pi/(1-nu)*uzA))
		uxDss=-((Dn0ss/4/pi/(1-nu)*uxA)+(Dss1__/4/pi/(1-nu)*uyA)+(Dds1ss/4/pi/(1-nu)*uzA))
		uxDds=-((Dn0ds/4/pi/(1-nu)*uxA)+(Dss0ds/4/pi/(1-nu)*uyA)+(Dds1__/4/pi/(1-nu)*uzA))
		uyDn =-((Dn1__/4/pi/(1-nu)*vxA)+(Dss0n_/4/pi/(1-nu)*vyA)+(Dds0n_/4/pi/(1-nu)*vzA))
		uyDss=-((Dn0ss/4/pi/(1-nu)*vxA)+(Dss1__/4/pi/(1-nu)*vyA)+(Dds1ss/4/pi/(1-nu)*vzA))
		uyDds=-((Dn0ds/4/pi/(1-nu)*vxA)+(Dss0ds/4/pi/(1-nu)*vyA)+(Dds1__/4/pi/(1-nu)*vzA))
		uzDn =-((Dn1__/4/pi/(1-nu)*wxA)+(Dss0n_/4/pi/(1-nu)*wyA)+(Dds0n_/4/pi/(1-nu)*wzA))
		uzDss=-((Dn0ss/4/pi/(1-nu)*wxA)+(Dss1__/4/pi/(1-nu)*wyA)+(Dds1ss/4/pi/(1-nu)*wzA))
		uzDds=-((Dn0ds/4/pi/(1-nu)*wxA)+(Dss0ds/4/pi/(1-nu)*wyA)+(Dds1__/4/pi/(1-nu)*wzA))
		
		
		(uxB,uyB,uzB,vxB,vyB,vzB,wxB,wyB,wzB) = AngDisDispFSC(y1B[indx[i]],y2B[indx[i]],y3B[indx[i]],cosB,sinB,cotB,cotB2,nu,-PB[3]);
		
		#Add mixed components together (different coords)
		uxDn  =uxDn  + ((Dn1__/4/pi/(1-nu)*uxB)+(Dss0n_/4/pi/(1-nu)*uyB)+(Dds0n_/4/pi/(1-nu)*uzB))
		uxDss =uxDss + ((Dn0ss/4/pi/(1-nu)*uxB)+(Dss1__/4/pi/(1-nu)*uyB)+(Dds1ss/4/pi/(1-nu)*uzB))
		uxDds =uxDds + ((Dn0ds/4/pi/(1-nu)*uxB)+(Dss0ds/4/pi/(1-nu)*uyB)+(Dds1__/4/pi/(1-nu)*uzB))
		uyDn  =uyDn  + ((Dn1__/4/pi/(1-nu)*vxB)+(Dss0n_/4/pi/(1-nu)*vyB)+(Dds0n_/4/pi/(1-nu)*vzB))
		uyDss =uyDss + ((Dn0ss/4/pi/(1-nu)*vxB)+(Dss1__/4/pi/(1-nu)*vyB)+(Dds1ss/4/pi/(1-nu)*vzB))
		uyDds =uyDds + ((Dn0ds/4/pi/(1-nu)*vxB)+(Dss0ds/4/pi/(1-nu)*vyB)+(Dds1__/4/pi/(1-nu)*vzB))
		uzDn  =uzDn  + ((Dn1__/4/pi/(1-nu)*wxB)+(Dss0n_/4/pi/(1-nu)*wyB)+(Dds0n_/4/pi/(1-nu)*wzB))
		uzDss =uzDss + ((Dn0ss/4/pi/(1-nu)*wxB)+(Dss1__/4/pi/(1-nu)*wyB)+(Dds1ss/4/pi/(1-nu)*wzB))
		uzDds =uzDds + ((Dn0ds/4/pi/(1-nu)*wxB)+(Dss0ds/4/pi/(1-nu)*wyB)+(Dds1__/4/pi/(1-nu)*wzB))
		
		#Add mixed components together (different coords)
		(uxDn, uyDn, uzDn) =RotateObject3DNewCoords!(uxDn, uyDn, uzDn ,0,0,0,VxR,VyR,VzR)
		(uxDss,uyDss,uzDss)=RotateObject3DNewCoords!(uxDss,uyDss,uzDss,0,0,0,VxR,VyR,VzR)
		(uxDds,uyDds,uzDds)=RotateObject3DNewCoords!(uxDds,uyDds,uzDds,0,0,0,VxR,VyR,VzR)
		
				
		#Add mixed components together (different coords)
		UxDn[indx[i]]  =UxDn[indx[i]]  + uxDn;
		UxDss[indx[i]] =UxDss[indx[i]] + uxDss;
		UxDds[indx[i]] =UxDds[indx[i]] + uxDds;
		UyDn[indx[i]]  =UyDn[indx[i]]  + uyDn;
		UyDss[indx[i]] =UyDss[indx[i]] + uyDss;
		UyDds[indx[i]] =UyDds[indx[i]] + uyDds;
		UzDn[indx[i]]  =UzDn[indx[i]]  + uzDn;
		UzDss[indx[i]] =UzDss[indx[i]] + uzDss;
		UzDds[indx[i]] =UzDds[indx[i]] + uzDds;
		
		
	end
	
	
	b=beta;
	sinB = sin(b);
	cosB = cos(b);
	cotB = cot(b);
	cotB2=cot(b/2);
	# Configuration II
	for i=1:length(indxf)
	
		(uxA,uyA,uzA,vxA,vyA,vzA,wxA,wyA,wzA) = AngDisDispFSC(y1A[indxf[i]],y2A[indxf[i]],y3A[indxf[i]],cosB,sinB,cotB,cotB2,nu,-PA[3]);
		
		#Add mixed components together (different coords)
		uxDn = -((Dn1__/4/pi/(1-nu)*uxA)+(Dss0n_/4/pi/(1-nu)*uyA)+(Dds0n_/4/pi/(1-nu)*uzA))
		uxDss= -((Dn0ss/4/pi/(1-nu)*uxA)+(Dss1__/4/pi/(1-nu)*uyA)+(Dds1ss/4/pi/(1-nu)*uzA))
		uxDds= -((Dn0ds/4/pi/(1-nu)*uxA)+(Dss0ds/4/pi/(1-nu)*uyA)+(Dds1__/4/pi/(1-nu)*uzA))
		uyDn = -((Dn1__/4/pi/(1-nu)*vxA)+(Dss0n_/4/pi/(1-nu)*vyA)+(Dds0n_/4/pi/(1-nu)*vzA))
		uyDss= -((Dn0ss/4/pi/(1-nu)*vxA)+(Dss1__/4/pi/(1-nu)*vyA)+(Dds1ss/4/pi/(1-nu)*vzA))
		uyDds= -((Dn0ds/4/pi/(1-nu)*vxA)+(Dss0ds/4/pi/(1-nu)*vyA)+(Dds1__/4/pi/(1-nu)*vzA))
		uzDn = -((Dn1__/4/pi/(1-nu)*wxA)+(Dss0n_/4/pi/(1-nu)*wyA)+(Dds0n_/4/pi/(1-nu)*wzA))
		uzDss= -((Dn0ss/4/pi/(1-nu)*wxA)+(Dss1__/4/pi/(1-nu)*wyA)+(Dds1ss/4/pi/(1-nu)*wzA))
		uzDds= -((Dn0ds/4/pi/(1-nu)*wxA)+(Dss0ds/4/pi/(1-nu)*wyA)+(Dds1__/4/pi/(1-nu)*wzA))	
		
		(uxB,uyB,uzB,vxB,vyB,vzB,wxB,wyB,wzB) = AngDisDispFSC(y1B[indxf[i]],y2B[indxf[i]],y3B[indxf[i]],cosB,sinB,cotB,cotB2,nu,-PB[3]);
		
		#Add mixed components together (different coords)
		uxDn  =uxDn  + ((Dn1__/4/pi/(1-nu)*uxB)+(Dss0n_/4/pi/(1-nu)*uyB)+(Dds0n_/4/pi/(1-nu)*uzB))
		uxDss =uxDss + ((Dn0ss/4/pi/(1-nu)*uxB)+(Dss1__/4/pi/(1-nu)*uyB)+(Dds1ss/4/pi/(1-nu)*uzB))
		uxDds =uxDds + ((Dn0ds/4/pi/(1-nu)*uxB)+(Dss0ds/4/pi/(1-nu)*uyB)+(Dds1__/4/pi/(1-nu)*uzB))
		uyDn  =uyDn  + ((Dn1__/4/pi/(1-nu)*vxB)+(Dss0n_/4/pi/(1-nu)*vyB)+(Dds0n_/4/pi/(1-nu)*vzB))
		uyDss =uyDss + ((Dn0ss/4/pi/(1-nu)*vxB)+(Dss1__/4/pi/(1-nu)*vyB)+(Dds1ss/4/pi/(1-nu)*vzB))
		uyDds =uyDds + ((Dn0ds/4/pi/(1-nu)*vxB)+(Dss0ds/4/pi/(1-nu)*vyB)+(Dds1__/4/pi/(1-nu)*vzB))
		uzDn  =uzDn  + ((Dn1__/4/pi/(1-nu)*wxB)+(Dss0n_/4/pi/(1-nu)*wyB)+(Dds0n_/4/pi/(1-nu)*wzB))
		uzDss =uzDss + ((Dn0ss/4/pi/(1-nu)*wxB)+(Dss1__/4/pi/(1-nu)*wyB)+(Dds1ss/4/pi/(1-nu)*wzB))
		uzDds =uzDds + ((Dn0ds/4/pi/(1-nu)*wxB)+(Dss0ds/4/pi/(1-nu)*wyB)+(Dds1__/4/pi/(1-nu)*wzB))
		
		#Add mixed components together (different coords)
		(uxDn, uyDn, uzDn) =RotateObject3DNewCoords!(uxDn, uyDn, uzDn ,0,0,0,VxR,VyR,VzR)
		(uxDss,uyDss,uzDss)=RotateObject3DNewCoords!(uxDss,uyDss,uzDss,0,0,0,VxR,VyR,VzR)
		(uxDds,uyDds,uzDds)=RotateObject3DNewCoords!(uxDds,uyDds,uzDds,0,0,0,VxR,VyR,VzR)
		
		#Add mixed components together (different coords)
		UxDn[indxf[i]]  =UxDn[indxf[i]]  + uxDn;
		UxDss[indxf[i]] =UxDss[indxf[i]] + uxDss;
		UxDds[indxf[i]] =UxDds[indxf[i]] + uxDds;
		UyDn[indxf[i]]  =UyDn[indxf[i]]  + uyDn;
		UyDss[indxf[i]] =UyDss[indxf[i]] + uyDss;
		UyDds[indxf[i]] =UyDds[indxf[i]] + uyDds;
		UzDn[indxf[i]]  =UzDn[indxf[i]]  + uzDn;
		UzDss[indxf[i]] =UzDss[indxf[i]] + uzDss;
		UzDds[indxf[i]] =UzDds[indxf[i]] + uzDds;
		
    end

end	#if statement

return( UxDn, UyDn, UzDn,
		UxDss,UyDss,UzDss,
		UxDds,UyDds,UzDds)
end


function CalcSideVec(PA,PB)
# Calculate TD side vector and the angle of the angular dislocation pair

SideVec = PB-PA;
eZ = [0;0;1];

G=-SideVec'*eZ/norm(SideVec);
beta = acos(G[1]);

return(SideVec,eZ,beta)
end

function CalcSlipVectorDiscCoords(SideVec,eZ,X,Y,Z,PA,beta)
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

# Determine the best arteact-free configuration for the calculation
# points near the free furface
I = (beta.*y1A).>=0;

return(I,y1A,y2A,y3A,y1B,y2B,y3B,ey1,ey2,ey3)
end


function AngDisDispFSC(y1,y2,y3,cosB,sinB,cotB,cotB2,nu,a)
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

#Export individual components
return(ux[1],uy[1],uz[1],vx[1],vy[1],vz[1],wx[1],wy[1],wz[1])

end


function TDstrainFS(X,Y,Z,P1,P2,P3,Dss,Dds,Dn,mu,lambda,ImageFlag,
		 ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
		 ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
		 ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds)
# TDstressFS 
# Calculates stresses and strains associated with a triangular dislocation 
# in an elastic full-space.

nu =lambda/(2*(mu+lambda));    #Poisson's ratio, Equation 8.28 Pollard

#Convert from global to the dislocation coordinate system
(Vnorm,Vstrike,Vdip)=CalculateLocalTriCoords(P1,P2,P3)
(p1,p2,p3,x,y,z,e12,e13,e23,A,B,C,Pos,Neg,casezLog)=GlobalToTDECoords(P1,P2,P3,X,Y,Z,Vnorm,Vstrike,Vdip)

if ImageFlag==1; #This means we are computing the iamge dislocation
	#Inverse rot mat
	VxR=[Vnorm[1],Vstrike[1],Vdip[1]];
	VyR=[Vnorm[2],Vstrike[2],Vdip[2]];
	VzR=[Vnorm[3],Vstrike[3],Vdip[3]];
	# Transform the strain tensor components from TDCS into EFCS
	(ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn) = TensorTransformation3D!(ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn,[VxR VyR VzR]);
	(ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss) = TensorTransformation3D!(ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,[VxR VyR VzR]);
	(ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds) = TensorTransformation3D!(ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,[VxR VyR VzR]);
	#We flip here and flip lower down (essentially taking away the components computed between these flips)
	if P1[3]==0 && P2[3]==0 && P3[3]==0
		ExzDn = -ExzDn;
		EyzDn = -EyzDn;
		ExzDss = -ExzDss;
		EyzDss = -EyzDss;
		ExzDds = -ExzDds;
		EyzDds = -EyzDds;
	end
end
	

#Turn to cart index
Pos=findall(Pos);
Neg=findall(Neg)

# Calculate first angular dislocation contribution POS
(ExxDn[Pos], EyyDn[Pos], EzzDn[Pos], ExyDn[Pos], ExzDn[Pos], EyzDn[Pos],
 ExxDss[Pos],EyyDss[Pos],EzzDss[Pos],ExyDss[Pos],ExzDss[Pos],EyzDss[Pos],
 ExxDds[Pos],EyyDds[Pos],EzzDds[Pos],ExyDds[Pos],ExzDds[Pos],EyzDds[Pos]) =
 TDSetupS(x[Pos],y[Pos],z[Pos],A,Dn,Dss,Dds,nu,p1,-e13,
 ExxDn[Pos], EyyDn[Pos], EzzDn[Pos], ExyDn[Pos], ExzDn[Pos], EyzDn[Pos],
 ExxDss[Pos],EyyDss[Pos],EzzDss[Pos],ExyDss[Pos],ExzDss[Pos],EyzDss[Pos],
 ExxDds[Pos],EyyDds[Pos],EzzDds[Pos],ExyDds[Pos],ExzDds[Pos],EyzDds[Pos]);
 
# Calculate second angular dislocation contribution
(ExxDn[Pos], EyyDn[Pos], EzzDn[Pos], ExyDn[Pos], ExzDn[Pos], EyzDn[Pos],
 ExxDss[Pos],EyyDss[Pos],EzzDss[Pos],ExyDss[Pos],ExzDss[Pos],EyzDss[Pos],
 ExxDds[Pos],EyyDds[Pos],EzzDds[Pos],ExyDds[Pos],ExzDds[Pos],EyzDds[Pos]) =
 TDSetupS(x[Pos],y[Pos],z[Pos],B,Dn,Dss,Dds,nu,p2,e12,
 ExxDn[Pos], EyyDn[Pos], EzzDn[Pos], ExyDn[Pos], ExzDn[Pos], EyzDn[Pos],
 ExxDss[Pos],EyyDss[Pos],EzzDss[Pos],ExyDss[Pos],ExzDss[Pos],EyzDss[Pos],
 ExxDds[Pos],EyyDds[Pos],EzzDds[Pos],ExyDds[Pos],ExzDds[Pos],EyzDds[Pos]); 
 
# Calculate third angular dislocation contribution
(ExxDn[Pos], EyyDn[Pos], EzzDn[Pos], ExyDn[Pos], ExzDn[Pos], EyzDn[Pos],
 ExxDss[Pos],EyyDss[Pos],EzzDss[Pos],ExyDss[Pos],ExzDss[Pos],EyzDss[Pos],
 ExxDds[Pos],EyyDds[Pos],EzzDds[Pos],ExyDds[Pos],ExzDds[Pos],EyzDds[Pos]) =
 TDSetupS(x[Pos],y[Pos],z[Pos],C,Dn,Dss,Dds,nu,p3,e23,
 ExxDn[Pos], EyyDn[Pos], EzzDn[Pos], ExyDn[Pos], ExzDn[Pos], EyzDn[Pos],
 ExxDss[Pos],EyyDss[Pos],EzzDss[Pos],ExyDss[Pos],ExzDss[Pos],EyzDss[Pos],
 ExxDds[Pos],EyyDds[Pos],EzzDds[Pos],ExyDds[Pos],ExzDds[Pos],EyzDds[Pos]);


# Calculate first angular dislocation contribution NEG
(ExxDn[Neg], EyyDn[Neg], EzzDn[Neg], ExyDn[Neg], ExzDn[Neg], EyzDn[Neg],
 ExxDss[Neg],EyyDss[Neg],EzzDss[Neg],ExyDss[Neg],ExzDss[Neg],EyzDss[Neg],
 ExxDds[Neg],EyyDds[Neg],EzzDds[Neg],ExyDds[Neg],ExzDds[Neg],EyzDds[Neg]) =
 TDSetupS(x[Neg],y[Neg],z[Neg],A,Dn,Dss,Dds,nu,p1,e13,
 ExxDn[Neg], EyyDn[Neg], EzzDn[Neg], ExyDn[Neg], ExzDn[Neg], EyzDn[Neg],
 ExxDss[Neg],EyyDss[Neg],EzzDss[Neg],ExyDss[Neg],ExzDss[Neg],EyzDss[Neg],
 ExxDds[Neg],EyyDds[Neg],EzzDds[Neg],ExyDds[Neg],ExzDds[Neg],EyzDds[Neg]);
 
# Calculate second angular dislocation contribution
(ExxDn[Neg], EyyDn[Neg], EzzDn[Neg], ExyDn[Neg], ExzDn[Neg], EyzDn[Neg],
 ExxDss[Neg],EyyDss[Neg],EzzDss[Neg],ExyDss[Neg],ExzDss[Neg],EyzDss[Neg],
 ExxDds[Neg],EyyDds[Neg],EzzDds[Neg],ExyDds[Neg],ExzDds[Neg],EyzDds[Neg]) = 
 TDSetupS(x[Neg],y[Neg],z[Neg],B,Dn,Dss,Dds,nu,p2,-e12,
 ExxDn[Neg], EyyDn[Neg], EzzDn[Neg], ExyDn[Neg], ExzDn[Neg], EyzDn[Neg],
 ExxDss[Neg],EyyDss[Neg],EzzDss[Neg],ExyDss[Neg],ExzDss[Neg],EyzDss[Neg],
 ExxDds[Neg],EyyDds[Neg],EzzDds[Neg],ExyDds[Neg],ExzDds[Neg],EyzDds[Neg]);
 
# Calculate third angular dislocation contribution 
(ExxDn[Neg], EyyDn[Neg], EzzDn[Neg], ExyDn[Neg], ExzDn[Neg], EyzDn[Neg],
 ExxDss[Neg],EyyDss[Neg],EzzDss[Neg],ExyDss[Neg],ExzDss[Neg],EyzDss[Neg],
 ExxDds[Neg],EyyDds[Neg],EzzDds[Neg],ExyDds[Neg],ExzDds[Neg],EyzDds[Neg]) = 
 TDSetupS(x[Neg],y[Neg],z[Neg],C,Dn,Dss,Dds,nu,p3,-e23,
 ExxDn[Neg], EyyDn[Neg], EzzDn[Neg], ExyDn[Neg], ExzDn[Neg], EyzDn[Neg],
 ExxDss[Neg],EyyDss[Neg],EzzDss[Neg],ExyDss[Neg],ExzDss[Neg],EyzDss[Neg],
 ExxDds[Neg],EyyDds[Neg],EzzDds[Neg],ExyDds[Neg],ExzDds[Neg],EyzDds[Neg]);	
		
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

 
# Transform the strain tensor components from TDCS into EFCS
(ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn) = TensorTransformation3D!(ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn,[Vnorm Vstrike Vdip]);
(ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss) = TensorTransformation3D!(ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,[Vnorm Vstrike Vdip]);
(ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds) = TensorTransformation3D!(ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,[Vnorm Vstrike Vdip]);

if ImageFlag==1; #This means we are computing the iamge dislocation
	if P1[3]==0 && P2[3]==0 && P3[3]==0
		ExzDn = -ExzDn;
		EyzDn = -EyzDn;
		ExzDss = -ExzDss;
		EyzDss = -EyzDss;
		ExzDds = -ExzDds;
		EyzDds = -EyzDds;
	end
end

return(	ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn,
		ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
		ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds)
 
end


function TDSetupS(x,y,z,alpha,Dn,Dss,Dds,nu,TriVertex,SideVec,
				 ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
				 ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
				 ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds)
# TDSetupS transforms coordinates of the calculation points as well as 
# slip vector components from ADCS into TDCS. It then calculates the 
# strains in ADCS and transforms them into TDCS.

Ct=SideVec[3];
St=SideVec[2];
# Transform coordinates of the calculation points from TDCS into ADCS
(y,z)  =RotateObject2D!(y,z,TriVertex[2],TriVertex[3],Ct,St)
# Transform the in-plane slip vector components from TDCS into ADCS
(Dss1,Dds0)=RotateObject2D(Dss,0.,0,0,Ct,St)
(Dss0,Dds1)=RotateObject2D(0.,Dds,0,0,Ct,St)

#Init arrays
Ang=-pi+alpha;
cosA = cos(Ang);
sinA = sin(Ang);

#Extra defs out of loop to speed it up
E1=(1-nu); #Elastic cons
E2=(2*nu+1);
cosA2=cosA^2;
sinADE1=sinA/8/pi/(1-nu);

#Burger func constants
E3=Dn/8/pi/E1;

E5Dss1=Dss1/8/pi/E1;
E5Dss0=Dss0/8/pi/E1;


# Transform strains from ADCS into TDCS 
B=[[1   0   0 ]; [0  Ct  St ]; [0  -St  Ct]]	 # 3x3 Transformation matrix

# Calculate strains associated with an angular dislocation in ADCS
for i=1:length(x)

	(exxDn,exxDss,exxDds,
	 eyyDn,eyyDss,eyyDds,
	 ezzDn,ezzDss,ezzDds,
	 exyDn,exyDss,exyDds,
	 exzDn,exzDss,exzDds,
	 eyzDn,eyzDss,eyzDds,
	 rFi_rx,rFi_ry,rFi_rz) =
	 AngDisStrain(x[i],y[i],z[i],cosA,sinA,Dn,Dss1[1],Dds1[1],nu,E1,E2,cosA2,sinADE1)

	#Opening components 
	exxdn= Dn*(rFi_rx)+E3* exxDn
	eyydn= E3* eyyDn; 
	ezzdn= E3* ezzDn;
	exydn= Dn*(rFi_ry)/2-E3* exyDn; 
	exzdn= Dn*(rFi_rz)/2-E3* exzDn;
	eyzdn= E3* eyzDn;
	
	#Some constants
	E4Dss1=Dss1*x[i]/8/pi/E1;
	E4Dss0=Dss0*x[i]/8/pi/E1;
	
	#Comp mixed components (in the current coords)
	exxdss= (-E4Dss1* exxDss)  +  (Dds0*x[i]*sinADE1* exxDds);
	exxdds= (-E4Dss0* exxDss)  +  (Dds1*x[i]*sinADE1* exxDds);
	
	eyydss= (Dss1*(rFi_ry)-E4Dss1* eyyDss)  +  (Dds0*x[i]*sinADE1* eyyDds);
	eyydds= (Dss0*(rFi_ry)-E4Dss0* eyyDss)  +  (Dds1*x[i]*sinADE1* eyyDds);
	
	ezzdss= (-E4Dss1* ezzDss) + (Dds0*(rFi_rz)+Dds0*x[i]*sinADE1* ezzDds); 
	ezzdds= (-E4Dss0* ezzDss) + (Dds1*(rFi_rz)+Dds1*x[i]*sinADE1* ezzDds); 
	
	exydss= (Dss1*(rFi_rx)/2+E5Dss1* exyDss)  -  (Dds0*sinADE1* exyDds);	
	exydds= (Dss0*(rFi_rx)/2+E5Dss0* exyDss)  -  (Dds1*sinADE1* exyDds);	
	
	exzdss= (E5Dss1* exzDss) + (Dds0*(rFi_rx)/2-Dds0*sinADE1*exzDds) ; 
	exzdds= (E5Dss0* exzDss) + (Dds1*(rFi_rx)/2-Dds1*sinADE1*exzDds);

	eyzdss= (Dss1*(rFi_rz)/2-E4Dss1* eyzDss) + (Dds0*(rFi_ry)/2-Dds0*x[i]*sinADE1* eyzDds);
	eyzdds= (Dss0*(rFi_rz)/2-E4Dss0* eyzDss) + (Dds1*(rFi_ry)/2-Dds1*x[i]*sinADE1* eyzDds);		  

	# Transform strains from ADCS into TDCS
	(exxdn,eyydn,ezzdn,exydn,exzdn,eyzdn) = 
	TensorTransformation3D!(exxdn,eyydn,ezzdn,exydn,exzdn,eyzdn,B);
	(exxdss,eyydss,ezzdss,exydss,exzdss,eyzdss) = 
	TensorTransformation3D!(exxdss,eyydss,ezzdss,exydss,exzdss,eyzdss,B);
	(exxdds,eyydds,ezzdds,exydds,exzdds,eyzdds) = 
	TensorTransformation3D!(exxdds,eyydds,ezzdds,exydds,exzdds,eyzdds,B);
	
	#Add to total vector
	ExxDn[i] 	= ExxDn[i] +exxdn;
	EyyDn[i] 	= EyyDn[i] +eyydn
	EzzDn[i] 	= EzzDn[i] +ezzdn
	ExyDn[i] 	= ExyDn[i] +exydn
	ExzDn[i] 	= ExzDn[i] +exzdn
	EyzDn[i] 	= EyzDn[i] +eyzdn
	
	ExxDss[i] 	= ExxDss[i] +exxdss;
	EyyDss[i] 	= EyyDss[i] +eyydss
	EzzDss[i] 	= EzzDss[i] +ezzdss
	ExyDss[i] 	= ExyDss[i] +exydss
	ExzDss[i] 	= ExzDss[i] +exzdss
	EyzDss[i] 	= EyzDss[i] +eyzdss
	
	ExxDds[i] 	= ExxDds[i] +exxdds
	EyyDds[i] 	= EyyDds[i] +eyydds
	EzzDds[i] 	= EzzDds[i] +ezzdds
	ExyDds[i] 	= ExyDds[i] +exydds
	ExzDds[i] 	= ExzDds[i] +exzdds
	EyzDds[i] 	= EyzDds[i] +eyzdds
	
	
end	


## Transform coordinates of the calculation points from ADCS into TDCS
(y,z)  =RotateObject2D!(y,z,0,0,Ct,-St)
y.=y.-TriVertex[2];
z.=z.-TriVertex[3];

return(	ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
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


function TDstrain_HarFunc(X,Y,Z,P1,P2,P3,Dss,Dds,Dn,mu,lambda,nu,
			 ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
			 ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
			 ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds);
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

# Calculate contribution of angular dislocation pair on each TD side 

# P1P2
(ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
 ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
 ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds) =
 AngSetupStrainFSC(X,Y,Z,Dn,Dss,Dds,P1,P2,mu,lambda,nu,Vnorm,Vstrike,Vdip,
 ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
 ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
 ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds);
 
# P2P3
(ExxDn13, EyyDn13, EzzDn13, ExyDn13, ExzDn13, EyzDn13,
 ExxDss13,EyyDss13,EzzDss13,ExyDss13,ExzDss13,EyzDss13,
 ExxDds13,EyyDds13,EzzDds13,ExyDds13,ExzDds13,EyzDds13) =
 AngSetupStrainFSC(X,Y,Z,Dn,Dss,Dds,P2,P3,mu,lambda,nu,Vnorm,Vstrike,Vdip,
 ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
 ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
 ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds); 
 
 
# P3P1 
(ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
 ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
 ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds) =
 AngSetupStrainFSC(X,Y,Z,Dn,Dss,Dds,P3,P1,mu,lambda,nu,Vnorm,Vstrike,Vdip,
 ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
 ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
 ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds);

return(	ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn,
		ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
		ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds);
end


function AngSetupStrainFSC(X,Y,Z,Dn,Dss,Dds,PA,PB,mu,lambda,nu,Vnorm,Vstrike,Vdip,
							ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
							ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
							ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds);
# AngSetupFSC_S calculates the Free Surface Correction to strains and 
# stresses associated with angular dislocation pair on each TD side.

# Calculate TD side vector and the angle of the angular dislocation pair
(SideVec,eZ,beta)=CalcSideVec(PA,PB)

if abs(beta)<eps() || abs(pi-beta)<eps()
    #Inputs come out the same (we add 0...)
else
    (I,y1A,y2A,y3A,y1B,y2B,y3B,ey1,ey2,ey3)=CalcSlipVectorDiscCoords(SideVec,eZ,X,Y,Z,PA,beta)

	AFlip = [ey1[1] ey1[2] ey1[3] ey2[1] ey2[2]  ey2[3]  ey3[1] ey3[2]  ey3[3]]; # Transformation matrix
	
	## Transform slip vector components from TDCS into EFCS
	(Dn1__,Dss0n_,Dds0n_) = RotateObject3DNewCoords(Dn,0.,0.,0,0,0,Vnorm,Vstrike,Vdip);
	(Dn0ss,Dss1__,Dds1ss) = RotateObject3DNewCoords(0.,Dss,0.,0,0,0,Vnorm,Vstrike,Vdip);
	(Dn0ds,Dss0ds,Dds1__) = RotateObject3DNewCoords(0.,0.,Dds,0,0,0,Vnorm,Vstrike,Vdip);
	# Transform slip vector components from EFCS to ADCS
	(Dn1__,Dss0n_,Dds0n_)=RotateObject3DNewCoords(Dn1__,Dss0n_,Dds0n_,0,0,0,ey1,ey2,ey3)
	(Dn0ss,Dss1__,Dds1ss)=RotateObject3DNewCoords(Dn0ss,Dss1__,Dds1ss,0,0,0,ey1,ey2,ey3)
	(Dn0ds,Dss0ds,Dds1__)=RotateObject3DNewCoords(Dn0ds,Dss0ds,Dds1__,0,0,0,ey1,ey2,ey3)
	
	
    Iflp=.!I; #Invert the bool
	indx=findall(I);
	indxf=findall(Iflp);
    
	#Init some vars
	b=pi-beta;
	sinB = sin(b);
	cosB = cos(b);
	cotB = cot(b);
	
    # Configuration I
	for i=1:length(indx)
		 #xx    yy    zz    xy    xz    yz
		(v11Ax,v22Ax,v33Ax,v12Ax,v13Ax,v23Ax,  #Dn
	     v11Ay,v22Ay,v33Ay,v12Ay,v13Ay,v23Ay,  #Dss
	     v11Az,v22Az,v33Az,v12Az,v13Az,v23Az) = AngDisStrainFSC(-y1A[indx[i]],-y2A[indx[i]],y3A[indx[i]],cosB,sinB,cotB,nu,-PA[3]);
		 
		v13Ax = -v13Ax;
		v13Ay = -v13Ay;
		v13Az = -v13Az;
		v23Ax = -v23Ax;
		v23Ay = -v23Ay;
		v23Az = -v23Az;
		
		#For configuration 1 Dn and Dss are flipped
		exxDn = -((-Dn1__*v11Ax)+(-Dss0n_*v11Ay)+(Dds0n_*v11Az))
		exxDss= -((-Dn0ss*v11Ax)+(-Dss1__*v11Ay)+(Dds1ss*v11Az))
		exxDds= -((-Dn0ds*v11Ax)+(-Dss0ds*v11Ay)+(Dds1__*v11Az))
		eyyDn = -((-Dn1__*v22Ax)+(-Dss0n_*v22Ay)+(Dds0n_*v22Az))
		eyyDss= -((-Dn0ss*v22Ax)+(-Dss1__*v22Ay)+(Dds1ss*v22Az))
		eyyDds= -((-Dn0ds*v22Ax)+(-Dss0ds*v22Ay)+(Dds1__*v22Az))
		ezzDn = -((-Dn1__*v33Ax)+(-Dss0n_*v33Ay)+(Dds0n_*v33Az))
		ezzDss= -((-Dn0ss*v33Ax)+(-Dss1__*v33Ay)+(Dds1ss*v33Az))
		ezzDds= -((-Dn0ds*v33Ax)+(-Dss0ds*v33Ay)+(Dds1__*v33Az))
		exyDn = -((-Dn1__/2*v12Ax)+(-Dss0n_/2*v12Ay)+(Dds0n_/2*v12Az))
		exyDss= -((-Dn0ss/2*v12Ax)+(-Dss1__/2*v12Ay)+(Dds1ss/2*v12Az))
		exyDds= -((-Dn0ds/2*v12Ax)+(-Dss0ds/2*v12Ay)+(Dds1__/2*v12Az))
		exzDn = -((-Dn1__/2*v13Ax)+(-Dss0n_/2*v13Ay)+(Dds0n_/2*v13Az))
		exzDss= -((-Dn0ss/2*v13Ax)+(-Dss1__/2*v13Ay)+(Dds1ss/2*v13Az))
		exzDds= -((-Dn0ds/2*v13Ax)+(-Dss0ds/2*v13Ay)+(Dds1__/2*v13Az))
		eyzDn = -((-Dn1__/2*v23Ax)+(-Dss0n_/2*v23Ay)+(Dds0n_/2*v23Az))
		eyzDss= -((-Dn0ss/2*v23Ax)+(-Dss1__/2*v23Ay)+(Dds1ss/2*v23Az))
		eyzDds= -((-Dn0ds/2*v23Ax)+(-Dss0ds/2*v23Ay)+(Dds1__/2*v23Az))
		
		(v11Bx,v22Bx,v33Bx,v12Bx,v13Bx,v23Bx,
	     v11By,v22By,v33By,v12By,v13By,v23By,
	     v11Bz,v22Bz,v33Bz,v12Bz,v13Bz,v23Bz) = AngDisStrainFSC(-y1B[indx[i]],-y2B[indx[i]],y3B[indx[i]],cosB,sinB,cotB,nu,-PB[3]);
		 
		 
		v13Bx = -v13Bx;
		v13By = -v13By;
		v13Bz = -v13Bz;
		v23Bx = -v23Bx;
		v23By = -v23By;
		v23Bz = -v23Bz;
		
		exxDn = exxDn + ((-Dn1__*v11Bx)+(-Dss0n_*v11By)+(Dds0n_*v11Bz))
		exxDss= exxDss+ ((-Dn0ss*v11Bx)+(-Dss1__*v11By)+(Dds1ss*v11Bz))
		exxDds= exxDds+ ((-Dn0ds*v11Bx)+(-Dss0ds*v11By)+(Dds1__*v11Bz))
		eyyDn = eyyDn + ((-Dn1__*v22Bx)+(-Dss0n_*v22By)+(Dds0n_*v22Bz))
		eyyDss= eyyDss+ ((-Dn0ss*v22Bx)+(-Dss1__*v22By)+(Dds1ss*v22Bz))
		eyyDds= eyyDds+ ((-Dn0ds*v22Bx)+(-Dss0ds*v22By)+(Dds1__*v22Bz))
		ezzDn = ezzDn + ((-Dn1__*v33Bx)+(-Dss0n_*v33By)+(Dds0n_*v33Bz))
		ezzDss= ezzDss+ ((-Dn0ss*v33Bx)+(-Dss1__*v33By)+(Dds1ss*v33Bz))
		ezzDds= ezzDds+ ((-Dn0ds*v33Bx)+(-Dss0ds*v33By)+(Dds1__*v33Bz))
		exyDn = exyDn + ((-Dn1__/2*v12Bx)+(-Dss0n_/2*v12By)+(Dds0n_/2*v12Bz))
		exyDss= exyDss+ ((-Dn0ss/2*v12Bx)+(-Dss1__/2*v12By)+(Dds1ss/2*v12Bz))
		exyDds= exyDds+ ((-Dn0ds/2*v12Bx)+(-Dss0ds/2*v12By)+(Dds1__/2*v12Bz))
		exzDn = exzDn + ((-Dn1__/2*v13Bx)+(-Dss0n_/2*v13By)+(Dds0n_/2*v13Bz))
		exzDss= exzDss+ ((-Dn0ss/2*v13Bx)+(-Dss1__/2*v13By)+(Dds1ss/2*v13Bz))
		exzDds= exzDds+ ((-Dn0ds/2*v13Bx)+(-Dss0ds/2*v13By)+(Dds1__/2*v13Bz))
		eyzDn = eyzDn + ((-Dn1__/2*v23Bx)+(-Dss0n_/2*v23By)+(Dds0n_/2*v23Bz))
		eyzDss= eyzDss+ ((-Dn0ss/2*v23Bx)+(-Dss1__/2*v23By)+(Dds1ss/2*v23Bz))
		eyzDds= eyzDds+ ((-Dn0ds/2*v23Bx)+(-Dss0ds/2*v23By)+(Dds1__/2*v23Bz))
		
		# Transform total Free Surface Correction to strains from ADCS to EFCS
		(exxDn,eyyDn,ezzDn,exyDn,exzDn,eyzDn) = TensorTransformation3D!(exxDn,eyyDn,ezzDn,exyDn,exzDn,eyzDn,AFlip);
		(exxDss,eyyDss,ezzDss,exyDss,exzDss,eyzDss) = TensorTransformation3D!(exxDss,eyyDss,ezzDss,exyDss,exzDss,eyzDss,AFlip);
		(exxDds,eyyDds,ezzDds,exyDds,exzDds,eyzDds) = TensorTransformation3D!(exxDds,eyyDds,ezzDds,exyDds,exzDds,eyzDds,AFlip);	
		
		#Add to total vector
		ExxDn[indx[i]] = ExxDn[indx[i]] +exxDn;
		ExxDss[indx[i]]= ExxDss[indx[i]]+exxDss;
		ExxDds[indx[i]]= ExxDds[indx[i]]+exxDds;
		EyyDn[indx[i]] = EyyDn[indx[i]] +eyyDn;
		EyyDss[indx[i]]= EyyDss[indx[i]]+eyyDss;
		EyyDds[indx[i]]= EyyDds[indx[i]]+eyyDds;
		EzzDn[indx[i]] = EzzDn[indx[i]] +ezzDn;
		EzzDss[indx[i]]= EzzDss[indx[i]]+ezzDss;
		EzzDds[indx[i]]= EzzDds[indx[i]]+ezzDds;
		ExyDn[indx[i]] = ExyDn[indx[i]] +exyDn;
		ExyDss[indx[i]]= ExyDss[indx[i]]+exyDss;
		ExyDds[indx[i]]= ExyDds[indx[i]]+exyDds;
		ExzDn[indx[i]] = ExzDn[indx[i]] +exzDn;
		ExzDss[indx[i]]= ExzDss[indx[i]]+exzDss;
		ExzDds[indx[i]]= ExzDds[indx[i]]+exzDds;
		EyzDn[indx[i]] = EyzDn[indx[i]] +eyzDn;
		EyzDss[indx[i]]= EyzDss[indx[i]]+eyzDss;
		EyzDds[indx[i]]= EyzDds[indx[i]]+eyzDds;
		
		
	end
    
	b=beta;
	sinB = sin(b);
	cosB = cos(b);
	cotB = cot(b);
    # Configuration II
	for i=1:length(indxf)
	

	
		(v11Ax,v22Ax,v33Ax,v12Ax,v13Ax,v23Ax,
	     v11Ay,v22Ay,v33Ay,v12Ay,v13Ay,v23Ay,
	     v11Az,v22Az,v33Az,v12Az,v13Az,v23Az) = AngDisStrainFSC(y1A[indxf[i]],y2A[indxf[i]],y3A[indxf[i]],cosB,sinB,cotB,nu,-PA[3]);
		 
		#For configuration 2 no components are flipped
		exxDn =-((Dn1__*v11Ax)+(Dss0n_*v11Ay)+(Dds0n_*v11Az))
		exxDss=-((Dn0ss*v11Ax)+(Dss1__*v11Ay)+(Dds1ss*v11Az))
		exxDds=-((Dn0ds*v11Ax)+(Dss0ds*v11Ay)+(Dds1__*v11Az))
		eyyDn =-((Dn1__*v22Ax)+(Dss0n_*v22Ay)+(Dds0n_*v22Az))
		eyyDss=-((Dn0ss*v22Ax)+(Dss1__*v22Ay)+(Dds1ss*v22Az))
		eyyDds=-((Dn0ds*v22Ax)+(Dss0ds*v22Ay)+(Dds1__*v22Az))
		ezzDn =-((Dn1__*v33Ax)+(Dss0n_*v33Ay)+(Dds0n_*v33Az))
		ezzDss=-((Dn0ss*v33Ax)+(Dss1__*v33Ay)+(Dds1ss*v33Az))
		ezzDds=-((Dn0ds*v33Ax)+(Dss0ds*v33Ay)+(Dds1__*v33Az))
		exyDn =-((Dn1__/2*v12Ax)+(Dss0n_/2*v12Ay)+(Dds0n_/2*v12Az))
		exyDss=-((Dn0ss/2*v12Ax)+(Dss1__/2*v12Ay)+(Dds1ss/2*v12Az))
		exyDds=-((Dn0ds/2*v12Ax)+(Dss0ds/2*v12Ay)+(Dds1__/2*v12Az))
		exzDn =-((Dn1__/2*v13Ax)+(Dss0n_/2*v13Ay)+(Dds0n_/2*v13Az))
		exzDss=-((Dn0ss/2*v13Ax)+(Dss1__/2*v13Ay)+(Dds1ss/2*v13Az))
		exzDds=-((Dn0ds/2*v13Ax)+(Dss0ds/2*v13Ay)+(Dds1__/2*v13Az))
		eyzDn =-((Dn1__/2*v23Ax)+(Dss0n_/2*v23Ay)+(Dds0n_/2*v23Az))
		eyzDss=-((Dn0ss/2*v23Ax)+(Dss1__/2*v23Ay)+(Dds1ss/2*v23Az))
		eyzDds=-((Dn0ds/2*v23Ax)+(Dss0ds/2*v23Ay)+(Dds1__/2*v23Az))
		
		(v11Bx,v22Bx,v33Bx,v12Bx,v13Bx,v23Bx,
	     v11By,v22By,v33By,v12By,v13By,v23By,
	     v11Bz,v22Bz,v33Bz,v12Bz,v13Bz,v23Bz) = AngDisStrainFSC(y1B[indxf[i]],y2B[indxf[i]],y3B[indxf[i]],cosB,sinB,cotB,nu,-PB[3]);
		 
		exxDn = exxDn + ((Dn1__*v11Bx)+(Dss0n_*v11By)+(Dds0n_*v11Bz))
		exxDss= exxDss+ ((Dn0ss*v11Bx)+(Dss1__*v11By)+(Dds1ss*v11Bz))
		exxDds= exxDds+ ((Dn0ds*v11Bx)+(Dss0ds*v11By)+(Dds1__*v11Bz))
		eyyDn = eyyDn + ((Dn1__*v22Bx)+(Dss0n_*v22By)+(Dds0n_*v22Bz))
		eyyDss= eyyDss+ ((Dn0ss*v22Bx)+(Dss1__*v22By)+(Dds1ss*v22Bz))
		eyyDds= eyyDds+ ((Dn0ds*v22Bx)+(Dss0ds*v22By)+(Dds1__*v22Bz))
		ezzDn = ezzDn + ((Dn1__*v33Bx)+(Dss0n_*v33By)+(Dds0n_*v33Bz))
		ezzDss= ezzDss+ ((Dn0ss*v33Bx)+(Dss1__*v33By)+(Dds1ss*v33Bz))
		ezzDds= ezzDds+ ((Dn0ds*v33Bx)+(Dss0ds*v33By)+(Dds1__*v33Bz))
		exyDn = exyDn + ((Dn1__/2*v12Bx)+(Dss0n_/2*v12By)+(Dds0n_/2*v12Bz))
		exyDss= exyDss+ ((Dn0ss/2*v12Bx)+(Dss1__/2*v12By)+(Dds1ss/2*v12Bz))
		exyDds= exyDds+ ((Dn0ds/2*v12Bx)+(Dss0ds/2*v12By)+(Dds1__/2*v12Bz))
		exzDn = exzDn + ((Dn1__/2*v13Bx)+(Dss0n_/2*v13By)+(Dds0n_/2*v13Bz))
		exzDss= exzDss+ ((Dn0ss/2*v13Bx)+(Dss1__/2*v13By)+(Dds1ss/2*v13Bz))
		exzDds= exzDds+ ((Dn0ds/2*v13Bx)+(Dss0ds/2*v13By)+(Dds1__/2*v13Bz))
		eyzDn = eyzDn + ((Dn1__/2*v23Bx)+(Dss0n_/2*v23By)+(Dds0n_/2*v23Bz))
		eyzDss= eyzDss+ ((Dn0ss/2*v23Bx)+(Dss1__/2*v23By)+(Dds1ss/2*v23Bz))
		eyzDds= eyzDds+ ((Dn0ds/2*v23Bx)+(Dss0ds/2*v23By)+(Dds1__/2*v23Bz))
		
	     # Transform total Free Surface Correction to strains from ADCS to EFCS
		(exxDn,eyyDn,ezzDn,exyDn,exzDn,eyzDn) = TensorTransformation3D!(exxDn,eyyDn,ezzDn,exyDn,exzDn,eyzDn,AFlip);
		(exxDss,eyyDss,ezzDss,exyDss,exzDss,eyzDss) = TensorTransformation3D!(exxDss,eyyDss,ezzDss,exyDss,exzDss,eyzDss,AFlip);
		(exxDds,eyyDds,ezzDds,exyDds,exzDds,eyzDds) = TensorTransformation3D!(exxDds,eyyDds,ezzDds,exyDds,exzDds,eyzDds,AFlip);	
		
				
		#Add to total vector
		ExxDn[indxf[i]] = ExxDn[indxf[i]] +exxDn;
		ExxDss[indxf[i]]= ExxDss[indxf[i]]+exxDss;
		ExxDds[indxf[i]]= ExxDds[indxf[i]]+exxDds;
		EyyDn[indxf[i]] = EyyDn[indxf[i]] +eyyDn;
		EyyDss[indxf[i]]= EyyDss[indxf[i]]+eyyDss;
		EyyDds[indxf[i]]= EyyDds[indxf[i]]+eyyDds;
		EzzDn[indxf[i]] = EzzDn[indxf[i]] +ezzDn;
		EzzDss[indxf[i]]= EzzDss[indxf[i]]+ezzDss;
		EzzDds[indxf[i]]= EzzDds[indxf[i]]+ezzDds;
		ExyDn[indxf[i]] = ExyDn[indxf[i]] +exyDn;
		ExyDss[indxf[i]]= ExyDss[indxf[i]]+exyDss;
		ExyDds[indxf[i]]= ExyDds[indxf[i]]+exyDds;
		ExzDn[indxf[i]] = ExzDn[indxf[i]] +exzDn;
		ExzDss[indxf[i]]= ExzDss[indxf[i]]+exzDss;
		ExzDds[indxf[i]]= ExzDds[indxf[i]]+exzDds;
		EyzDn[indxf[i]] = EyzDn[indxf[i]] +eyzDn;
		EyzDss[indxf[i]]= EyzDss[indxf[i]]+eyzDss;
		EyzDds[indxf[i]]= EyzDds[indxf[i]]+eyzDds;
		
		
	end


	

end

return(	ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn,
		ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
		ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds);
end


function AngDisStrainFSC(y1,y2,y3,cosB,sinB,cotB,nu,a)
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

v11x = (1/4*((-2+2*nu)*N1*rFib_ry1*cotB^2-N1*y2/W6^2*((1-W5)*cotB-
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
    rb2^2*y1))/pi/(1-nu));
v11y = (1/4*(N1*(((2-2*nu)*cotB^2+nu)/rb*y1/W6-((2-2*nu)*cotB^2+1)*
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
    W9*(y1/rb-sinB)+W1/W7*a/rb^3/cosB*y1)))/pi/(1-nu));
v11z= (1/4*(N1*(-y2/W6^2*(1+a/rb)/rb*y1-y2/W6*a/rb^3*y1+y2*
    cosB/W7^2*W2*(y1/rb-sinB)+y2*cosB/W7*a/rb^3*y1)+y2*W8/
    rb^3*(a/rb2+1/W6)*y1-y2*W8/rb*(-2*a/rb2^2*y1-1/W6^2/
    rb*y1)-y2*W8*cosB/rb^3/W7*(W1/W7*W2+a*y3b/rb2)*y1-y2*W8*
    cosB/rb/W7^2*(W1/W7*W2+a*y3b/rb2)*(y1/rb-sinB)+y2*W8*
    cosB/rb/W7*(1/rb*cosB*y1/W7*W2-W1/W7^2*W2*(y1/rb-sinB)-
    W1/W7*a/rb^3*y1-2*a*y3b/rb2^2*y1))/pi/(1-nu));
	
v22x = (1/4*(N1*(((2-2*nu)*cotB^2-nu)/rb*y2/W6-((2-2*nu)*cotB^2+1-
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
    z1b*cotB/rb2*cosB*y2)))/pi/(1-nu));
v22y = (1/4*((2-2*nu)*N1*rFib_ry2*cotB^2+N1/W6*((W5-1)*cotB+y1/W6*
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
    2*a*y3b/rb2^2/cosB*y2))/pi/(1-nu));
v22z = (1/4*(N1*(-sinB/rb*y2/W7+y2/W6^2*(1+a/rb)/rb*y1+y2/W6*
    a/rb^3*y1-z1b/W7^2*W2/rb*y2-z1b/W7*a/rb^3*y2)-y2*W8/
    rb^3*(a/rb2+1/W6)*y1+y1*W8/rb*(-2*a/rb2^2*y2-1/W6^2/
    rb*y2)+W8/W7^2*(sinB*(cosB-a/rb)+z1b/rb*(1+a*y3b/rb2)-1/
    rb/W7*(y2^2*cosB*sinB-a*z1b/rb*W1))/rb*y2-W8/W7*(sinB*a/
    rb^3*y2-z1b/rb^3*(1+a*y3b/rb2)*y2-2*z1b/rb^5*a*y3b*y2+
    1/rb^3/W7*(y2^2*cosB*sinB-a*z1b/rb*W1)*y2+1/rb2/W7^2*
    (y2^2*cosB*sinB-a*z1b/rb*W1)*y2-1/rb/W7*(2*y2*cosB*sinB+a*
    z1b/rb^3*W1*y2-a*z1b/rb2*cosB*y2)))/pi/(1-nu));

v33x = (1/4*((2-2*nu)*(N1*rFib_ry3*cotB-y2/W6^2*W5*(y3b/rb+1)-
    1/2*y2/W6*a/rb^3*2*y3b+y2*cosB/W7^2*W2*W3+1/2*y2*cosB/W7*
    a/rb^3*2*y3b)+y2/rb*(2*nu/W6+a/rb2)-1/2*y2*W8/rb^3*(2*
    nu/W6+a/rb2)*2*y3b+y2*W8/rb*(-2*nu/W6^2*(y3b/rb+1)-a/
    rb2^2*2*y3b)+y2*cosB/rb/W7*(1-2*nu-W1/W7*W2-a*y3b/rb2)-
    1/2*y2*W8*cosB/rb^3/W7*(1-2*nu-W1/W7*W2-a*y3b/rb2)*2*
    y3b-y2*W8*cosB/rb/W7^2*(1-2*nu-W1/W7*W2-a*y3b/rb2)*W3+y2*
    W8*cosB/rb/W7*(-(cosB*y3b/rb+1)/W7*W2+W1/W7^2*W2*W3+1/2*
    W1/W7*a/rb^3*2*y3b-a/rb2+a*y3b/rb2^2*2*y3b))/pi/(1-nu));
v33y = (1/4*((-2+2*nu)*N1*cotB*((y3b/rb+1)/W6-cosB*W3/W7)+(2-2*nu)*
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
    W7*2*y3b+z1b*W1/rb/W7^2*W3)))/pi/(1-nu));
v33z =(1/4*((2-2*nu)*rFib_ry3-(2-2*nu)*y2*sinB/W7^2*W2*W3-1/2*
    (2-2*nu)*y2*sinB/W7*a/rb^3*2*y3b+y2*sinB/rb/W7*(1+W1/W7*
    W2+a*y3b/rb2)-1/2*y2*W8*sinB/rb^3/W7*(1+W1/W7*W2+a*y3b/
    rb2)*2*y3b-y2*W8*sinB/rb/W7^2*(1+W1/W7*W2+a*y3b/rb2)*W3+
    y2*W8*sinB/rb/W7*((cosB*y3b/rb+1)/W7*W2-W1/W7^2*W2*W3-
    1/2*W1/W7*a/rb^3*2*y3b+a/rb2-a*y3b/rb2^2*2*y3b))/pi/(1-nu));

v12x = (1/4*((-2+2*nu)*N1*rFib_ry2*cotB^2+N1/W6*((1-W5)*cotB-y1/
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
	(1/4*(N1*(((2-2*nu)*cotB^2-nu)/rb*y1/W6-((2-2*nu)*cotB^2+1-   
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
    rb2*cosB*y1)))/pi/(1-nu));
	
v12y = (1/4*(N1*(((2-2*nu)*cotB^2+nu)/rb*y2/W6-((2-2*nu)*cotB^2+1)*
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
	(1/4*((2-2*nu)*N1*rFib_ry1*cotB^2-N1*y2/W6^2*((W5-1)*cotB+
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
    W1/W7*a/rb^3/cosB*y1-2*a*y3b/rb2^2/cosB*y1))/pi/(1-nu));
	
v12z = (1/4*(N1*(1/W6*(1+a/rb)-y2^2/W6^2*(1+a/rb)/rb-y2^2/
    W6*a/rb^3-cosB/W7*W2+y2^2*cosB/W7^2*W2/rb+y2^2*cosB/W7*
    a/rb^3)-W8/rb*(a/rb2+1/W6)+y2^2*W8/rb^3*(a/rb2+1/W6)-
    y2*W8/rb*(-2*a/rb2^2*y2-1/W6^2/rb*y2)+W8*cosB/rb/W7*
    (W1/W7*W2+a*y3b/rb2)-y2^2*W8*cosB/rb^3/W7*(W1/W7*W2+a*
    y3b/rb2)-y2^2*W8*cosB/rb2/W7^2*(W1/W7*W2+a*y3b/rb2)+y2*
    W8*cosB/rb/W7*(1/rb*cosB*y2/W7*W2-W1/W7^2*W2/rb*y2-W1/
    W7*a/rb^3*y2-2*a*y3b/rb2^2*y2))/pi/(1-nu))+  
    (1/4*(N1*(-sinB*(y1/rb-sinB)/W7-1/W6*(1+a/rb)+y1^2/W6^
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

v13x = (1/4*((-2+2*nu)*N1*rFib_ry3*cotB^2-N1*y2/W6^2*((1-W5)*
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
	(1/4*((2-2*nu)*(N1*rFib_ry1*cotB-y1/W6^2*W5/rb*y2-y2/W6*
    a/rb^3*y1+y2*cosB/W7^2*W2*(y1/rb-sinB)+y2*cosB/W7*a/rb^
    3*y1)-y2*W8/rb^3*(2*nu/W6+a/rb2)*y1+y2*W8/rb*(-2*nu/W6^
    2/rb*y1-2*a/rb2^2*y1)-y2*W8*cosB/rb^3/W7*(1-2*nu-W1/W7*
    W2-a*y3b/rb2)*y1-y2*W8*cosB/rb/W7^2*(1-2*nu-W1/W7*W2-a*
    y3b/rb2)*(y1/rb-sinB)+y2*W8*cosB/rb/W7*(-1/rb*cosB*y1/W7*
    W2+W1/W7^2*W2*(y1/rb-sinB)+W1/W7*a/rb^3*y1+2*a*y3b/rb2^
    2*y1))/pi/(1-nu));
	
v13y = (1/4*(N1*(((2-2*nu)*cotB^2+nu)*(y3b/rb+1)/W6-((2-2*nu)*cotB^
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
	(1/4*((-2+2*nu)*N1*cotB*(1/rb*y1/W6-cosB*(y1/rb-sinB)/W7)-
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
    W1/rb/W7^2*(y1/rb-sinB))))/pi/(1-nu));
	
v13z = (1/4*(N1*(-y2/W6^2*(1+a/rb)*(y3b/rb+1)-1/2*y2/W6*a/
    rb^3*2*y3b+y2*cosB/W7^2*W2*W3+1/2*y2*cosB/W7*a/rb^3*2*
    y3b)-y2/rb*(a/rb2+1/W6)+1/2*y2*W8/rb^3*(a/rb2+1/W6)*2*
    y3b-y2*W8/rb*(-a/rb2^2*2*y3b-1/W6^2*(y3b/rb+1))+y2*cosB/
    rb/W7*(W1/W7*W2+a*y3b/rb2)-1/2*y2*W8*cosB/rb^3/W7*(W1/
    W7*W2+a*y3b/rb2)*2*y3b-y2*W8*cosB/rb/W7^2*(W1/W7*W2+a*
    y3b/rb2)*W3+y2*W8*cosB/rb/W7*((cosB*y3b/rb+1)/W7*W2-W1/
    W7^2*W2*W3-1/2*W1/W7*a/rb^3*2*y3b+a/rb2-a*y3b/rb2^2*2*
    y3b))/pi/(1-nu))+
    (1/4*((2-2*nu)*rFib_ry1-(2-2*nu)*y2*sinB/W7^2*W2*(y1/rb-
    sinB)-(2-2*nu)*y2*sinB/W7*a/rb^3*y1-y2*W8*sinB/rb^3/W7*(1+
    W1/W7*W2+a*y3b/rb2)*y1-y2*W8*sinB/rb/W7^2*(1+W1/W7*W2+
    a*y3b/rb2)*(y1/rb-sinB)+y2*W8*sinB/rb/W7*(1/rb*cosB*y1/
    W7*W2-W1/W7^2*W2*(y1/rb-sinB)-W1/W7*a/rb^3*y1-2*a*y3b/
    rb2^2*y1))/pi/(1-nu));

v23x = (1/4*(N1*(((2-2*nu)*cotB^2-nu)*(y3b/rb+1)/W6-((2-2*nu)*
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
	(1/4*((2-2*nu)*(N1*rFib_ry2*cotB+1/W6*W5-y2^2/W6^2*W5/
    rb-y2^2/W6*a/rb^3-cosB/W7*W2+y2^2*cosB/W7^2*W2/rb+y2^2*
    cosB/W7*a/rb^3)+W8/rb*(2*nu/W6+a/rb2)-y2^2*W8/rb^3*(2*
    nu/W6+a/rb2)+y2*W8/rb*(-2*nu/W6^2/rb*y2-2*a/rb2^2*y2)+
    W8*cosB/rb/W7*(1-2*nu-W1/W7*W2-a*y3b/rb2)-y2^2*W8*cosB/
    rb^3/W7*(1-2*nu-W1/W7*W2-a*y3b/rb2)-y2^2*W8*cosB/rb2/W7^
    2*(1-2*nu-W1/W7*W2-a*y3b/rb2)+y2*W8*cosB/rb/W7*(-1/rb*
    cosB*y2/W7*W2+W1/W7^2*W2/rb*y2+W1/W7*a/rb^3*y2+2*a*
    y3b/rb2^2*y2))/pi/(1-nu));
	
v23y = (1/4*((2-2*nu)*N1*rFib_ry3*cotB^2-N1*y2/W6^2*((W5-1)*cotB+
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
	(1/4*((-2+2*nu)*N1*cotB*(1/rb*y2/W6-cosB/rb*y2/W7)+(2-
    2*nu)*y1/W6^2*W5/rb*y2+(2-2*nu)*y1/W6*a/rb^3*y2-(2-2*
    nu)*z1b/W7^2*W2/rb*y2-(2-2*nu)*z1b/W7*a/rb^3*y2-W8/rb^
    3*(N1*cotB-2*nu*y1/W6-a*y1/rb2)*y2+W8/rb*(2*nu*y1/W6^2/
    rb*y2+2*a*y1/rb2^2*y2)+W8/W7^2*(cosB*sinB+W1*cotB/rb*((2-
    2*nu)*cosB-W1/W7)+a/rb*(sinB-y3b*z1b/rb2-z1b*W1/rb/W7))/
    rb*y2-W8/W7*(1/rb2*cosB*y2*cotB*((2-2*nu)*cosB-W1/W7)-W1*
    cotB/rb^3*((2-2*nu)*cosB-W1/W7)*y2+W1*cotB/rb*(-cosB/rb*
    y2/W7+W1/W7^2/rb*y2)-a/rb^3*(sinB-y3b*z1b/rb2-z1b*W1/
    rb/W7)*y2+a/rb*(2*y3b*z1b/rb2^2*y2-z1b/rb2*cosB*y2/W7+
    z1b*W1/rb^3/W7*y2+z1b*W1/rb2/W7^2*y2)))/pi/(1-nu));
	
v23z = (1/4*(N1*(-sinB*W3/W7+y1/W6^2*(1+a/rb)*(y3b/rb+1)+
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
    (1/4*((2-2*nu)*rFib_ry2+(2-2*nu)*sinB/W7*W2-(2-2*nu)*y2^2*
    sinB/W7^2*W2/rb-(2-2*nu)*y2^2*sinB/W7*a/rb^3+W8*sinB/rb/
    W7*(1+W1/W7*W2+a*y3b/rb2)-y2^2*W8*sinB/rb^3/W7*(1+W1/
    W7*W2+a*y3b/rb2)-y2^2*W8*sinB/rb2/W7^2*(1+W1/W7*W2+a*
    y3b/rb2)+y2*W8*sinB/rb/W7*(1/rb*cosB*y2/W7*W2-W1/W7^2*
    W2/rb*y2-W1/W7*a/rb^3*y2-2*a*y3b/rb2^2*y2))/pi/(1-nu));
	

return(v11x[1],v22x[1],v33x[1],v12x[1],v13x[1],v23x[1],
	   v11y[1],v22y[1],v33y[1],v12y[1],v13y[1],v23y[1],
	   v11z[1],v22z[1],v33z[1],v12z[1],v13z[1],v23z[1])
end