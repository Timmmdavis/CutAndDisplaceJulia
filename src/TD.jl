function TD(X::Array{Float64,2},Y::Array{Float64,2},Z::Array{Float64,2},
		    P1List::Array{Float64,2},P2List::Array{Float64,2},P3List::Array{Float64,2},
			Dss::Array,Dds::Array,Dn::Array,nu::Float64,mu::Float64,
			DispFlag::Int64,StrainFlag::Int64,HSflag::Int64)
# TD 
# Julia version of Mehdi's triangular dislocation functions, creates influence matrices. 
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
# Reference jthenal article: 
# Nikkhoo M. and Walter T.R., 2015. Triangular dislocation: An analytical, 
# artefact-free solution. 
# Submitted to Geophysical Jthenal International 
#
# Copyright (c) 2014 Mehdi Nikkhoo
# Copyright (c) 2018 Tim Davis
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

#Allocate some arrays of zeros that will be taken through the following loops with the results appended
if DispFlag==1
	UxDn  = zeros(length(X),SzCmp); 
	UyDn  = zeros(length(X),SzCmp); 
	UzDn  = zeros(length(X),SzCmp); 
	UxDss = zeros(length(X),SzCmp); 
	UyDss = zeros(length(X),SzCmp); 
	UzDss = zeros(length(X),SzCmp); 
	UxDds = zeros(length(X),SzCmp); 
	UyDds = zeros(length(X),SzCmp); 
	UzDds = zeros(length(X),SzCmp); 
end
if StrainFlag==1
	ExxDn  = zeros(length(X),SzCmp); 
	EyyDn  = zeros(length(X),SzCmp); 
	EzzDn  = zeros(length(X),SzCmp);
	ExyDn  = zeros(length(X),SzCmp); 
	ExzDn  = zeros(length(X),SzCmp); 
	EyzDn  = zeros(length(X),SzCmp);
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

P1iList=deepcopy(P1List);
P2iList=deepcopy(P2List);
P3iList=deepcopy(P3List);
for i=1:SzCmp
	P1iList[i,3] = -P1iList[i,3]
	P2iList[i,3] = -P2iList[i,3]
	P3iList[i,3] = -P3iList[i,3]
end


#Comp FaceNormalVector...

VnormList  	= zeros(SzCmp,3); 
VstrikeList = zeros(SzCmp,3); 
VdipList  	= zeros(SzCmp,3); 
VnormiList  	= zeros(SzCmp,3); 
VstrikeiList = zeros(SzCmp,3); 
VdipiList  	= zeros(SzCmp,3); 

#Some allocations out of loop
eY = [0.;1.;0.];
eZ = [0.;0.;1.];
FillA=zeros(3);
FillB=zeros(3);


Threads.@threads for i=1:SzCmp #For every element (multithreaded) 
	#println("Multithreading off")

	P1=view(P1List,i,:);
	P2=view(P2List,i,:);
	P3=view(P3List,i,:);
	P1i=view(P1iList,i,:);
	P2i=view(P2iList,i,:);
	P3i=view(P3iList,i,:);

	Vnorm=	view(VnormList,i,:);
	Vstrike=view(VstrikeList,i,:);
	Vdip=	view(VdipList,i,:);
	Vnormi=	view(VnormiList,i,:);
	Vstrikei=view(VstrikeiList,i,:);
	Vdipi=	view(VdipiList,i,:);
	
	CalculateLocalTriCoords!(P1,P2,P3,Vnorm,Vstrike,Vdip,eY,eZ,FillA,FillB);
	CalculateLocalTriCoords!(P1i,P2i,P3i,Vnormi,Vstrikei,Vdipi,eY,eZ,FillA,FillB);
	#error("Stop here, also add back in threads and ! funcs no output. ")
	
	#Allocate outside of funcs (using @view we just assign a pointer). 
	#See Gotcha #5 https://www.juliabloggers.com/7-julia-gotchas-and-how-to-handle-them/
	if DispFlag==1
		UxDnI  = view(UxDn,:,i); #only 48bytes alloc, note its the same as @view UxDn[:,i];
		UyDnI  = view(UyDn,:,i);
		UzDnI  = view(UzDn,:,i);
		UxDssI = view(UxDss,:,i);
		UyDssI = view(UyDss,:,i);
		UzDssI = view(UzDss,:,i);
		UxDdsI = view(UxDds,:,i);
		UyDdsI = view(UyDds,:,i);
		UzDdsI = view(UzDds,:,i);
	end
	if StrainFlag==1
		ExxDnI  = view(ExxDn,:,i);
		EyyDnI  = view(EyyDn,:,i);
		EzzDnI  = view(EzzDn,:,i);
		ExyDnI  = view(ExyDn,:,i);
		ExzDnI  = view(ExzDn,:,i);
		EyzDnI  = view(EyzDn,:,i);
		ExxDssI = view(ExxDss,:,i);
		EyyDssI = view(EyyDss,:,i);
		EzzDssI = view(EzzDss,:,i);
		ExyDssI = view(ExyDss,:,i);
		ExzDssI = view(ExzDss,:,i);
		EyzDssI = view(EyzDss,:,i);
		ExxDdsI = view(ExxDds,:,i);
		EyyDdsI = view(EyyDds,:,i);
		EzzDdsI = view(EzzDds,:,i);
		ExyDdsI = view(ExyDds,:,i);
		ExzDdsI = view(ExzDds,:,i);
		EyzDdsI = view(EyzDds,:,i);
	end

	if HSflag==1
		if any(Z .>0) | any(P1[3] .>0) | any(P2[3] .>0) | any(P3[3] .>0)
			error("Half-space solution: Z coordinates must be negative!")
		end
	end
	
	#Compute some variables we use repeated times inside the functions.
	##Reducing Allocs further in these would be good!
	(p1,p2,p3,x,y,z,e12,e13,e23,A,B,C,casepLog,casenLog,casezLog)=
	GlobalToTDECoords(P1,P2,P3,X,Y,Z,Vnorm,Vstrike,Vdip)
	(p1i,p2i,p3i,xi,yi,zi,e12i,e13i,e23i,Ai,Bi,Ci,casepLogi,casenLogi,casezLogi)=
	GlobalToTDECoords(P1i,P2i,P3i,X,Y,Z,Vnormi,Vstrikei,Vdipi)

	if DispFlag==1

		# Calculate main dislocation contribution to displacements
		TDdispFS(X,Y,Z,P1,P2,P3,Dss[i],Dds[i],Dn[i],nu,0,Vnorm,Vstrike,Vdip,
		p1,p2,p3,x,y,z,e12,e13,e23,A,B,C,casepLog,casenLog,casezLog,
		UxDnI, UyDnI, UzDnI,
		UxDssI,UyDssI,UzDssI,
		UxDdsI,UyDdsI,UzDdsI);
		 
		if HSflag==1
			
			# Calculate harmonic function contribution to displacements
			TDdisp_HarFunc(X,Y,Z,P1,P2,P3,Dss[i],Dds[i],Dn[i],nu,Vnorm,Vstrike,Vdip,
			UxDnI, UyDnI, UzDnI,
			UxDssI,UyDssI,UzDssI,
			UxDdsI,UyDdsI,UzDdsI);

			# Calculate image dislocation contribution to displacements
			TDdispFS(X,Y,Z,P1i,P2i,P3i,Dss[i],Dds[i],Dn[i],nu,1,Vnormi,Vstrikei,Vdipi,
			p1i,p2i,p3i,xi,yi,zi,e12i,e13i,e23i,Ai,Bi,Ci,casepLogi,casenLogi,casezLogi,
			UxDnI, UyDnI, UzDnI,
			UxDssI,UyDssI,UzDssI,
			UxDdsI,UyDdsI,UzDdsI);
	
		end #HsFlag

	end#DispFlag 

	if StrainFlag==1

		#Elastic con
		lambda=(2*mu*nu)/(1-(2*nu));
		
		# Calculate main dislocation contribution to strains
		TDstrainFS(X,Y,Z,P1,P2,P3,Dss[i],Dds[i],Dn[i],mu,lambda,nu,0,Vnorm,Vstrike,Vdip,
		p1,p2,p3,x,y,z,e12,e13,e23,A,B,C,casepLog,casenLog,casezLog,
		ExxDnI, EyyDnI, EzzDnI, ExyDnI, ExzDnI, EyzDnI,
		ExxDssI,EyyDssI,EzzDssI,ExyDssI,ExzDssI,EyzDssI,
		ExxDdsI,EyyDdsI,EzzDdsI,ExyDdsI,ExzDdsI,EyzDdsI);	
		
		if HSflag==1
		
			# Calculate harmonic function contribution to strains
			TDstrain_HarFunc(X,Y,Z,P1,P2,P3,Dss[i],Dds[i],Dn[i],mu,lambda,nu,Vnorm,Vstrike,Vdip,
			ExxDnI, EyyDnI, EzzDnI, ExyDnI, ExzDnI, EyzDnI,
			ExxDssI,EyyDssI,EzzDssI,ExyDssI,ExzDssI,EyzDssI,
			ExxDdsI,EyyDdsI,EzzDdsI,ExyDdsI,ExzDdsI,EyzDdsI);	
		
			# Calculate image dislocation contribution to strains 
			TDstrainFS(X,Y,Z,P1i,P2i,P3i,Dss[i],Dds[i],Dn[i],mu,lambda,nu,1,Vnormi,Vstrikei,Vdipi,
			p1i,p2i,p3i,xi,yi,zi,e12i,e13i,e23i,Ai,Bi,Ci,casepLogi,casenLogi,casezLogi,
			ExxDnI, EyyDnI, EzzDnI, ExyDnI, ExzDnI, EyzDnI,
			ExxDssI,EyyDssI,EzzDssI,ExyDssI,ExzDssI,EyzDssI,
			ExxDdsI,EyyDdsI,EzzDdsI,ExyDdsI,ExzDdsI,EyzDdsI);


		end #HS flag
		
	end #Stress flag
	
end #Over all els
	
	
if DispFlag==0
	UxDn  = [];
	UyDn  = [];
	UzDn  = [];
	UxDss = [];
	UyDss = [];
	UzDss = [];
	UxDds = [];
	UyDds = [];
	UzDds = [];
end
if StrainFlag==0
	ExxDn  = []; 
	EyyDn  = []; 
	EzzDn  = [];
	ExyDn  = []; 
	ExzDn  = []; 
	EyzDn  = [];
	ExxDss = []; 
	EyyDss = []; 
	EzzDss = [];
	ExyDss = []; 
	ExzDss = []; 
	EyzDss = [];
	ExxDds = []; 
	EyyDds = []; 
	EzzDds = [];
	ExyDds = []; 
	ExzDds = []; 
	EyzDds = [];
end

return(ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
	   ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
	   ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,
	   UxDn,UyDn,UzDn,
	   UxDss,UyDss,UzDss,
	   UxDds,UyDds,UzDds)
	

end

function CalculateLocalTriCoords!(P1,P2,P3,Vnorm,Vstrike,Vdip,eY,eZ,P1P2,P3P1)
# Calculate unit strike, dip and normal to TD vectors: For a horizontal TD 
# as an exception, if the normal vector points upward, the strike and dip 
# vectors point Northward and Westward, whereas if the normal vector points
# downward, the strike and dip vectors point Southward and Westward, 
# respectively.

P1P2.=P2.-P1;
P3P1.=P3.-P1;
cross!(P1P2,P3P1,Vnorm);
Vnorm.= Vnorm./sqrt(Vnorm[1]^2+Vnorm[2]^2+Vnorm[3]^2)
cross!(eZ,Vnorm,Vstrike);
if norm(Vstrike)==0
	Vstrike.= eY.*Vnorm[3];
end
Vstrike.= Vstrike./sqrt(Vstrike[1]^2+Vstrike[2]^2+Vstrike[3]^2);
cross!(Vnorm,Vstrike,Vdip);

end


function TDdispFS(X,Y,Z,P1,P2,P3,Dss,Dds,Dn,nu,ImageFlag,Vnorm,Vstrike,Vdip,
				p1,p2,p3,x,y,z,e12,e13,e23,A,B,C,Pos,Neg,casezLog,
				UxDn, UyDn, UzDn,
				UxDss,UyDss,UzDss,
				UxDds,UyDds,UzDds)


if ImageFlag==1; #This means we are computing the iamge dislocation

	#Inverse rot mat
	VxR=zeros(size(Vnorm));#Single alloc here!
	VyR=zeros(size(Vnorm));
	VzR=zeros(size(Vnorm));
	VxR[1]=Vnorm[1]; VxR[2]=Vstrike[1]; VxR[3]=Vdip[1];
	VyR[1]=Vnorm[2]; VyR[2]=Vstrike[2]; VyR[3]=Vdip[2];
	VzR[1]=Vnorm[3]; VzR[2]=Vstrike[3]; VzR[3]=Vdip[3];

	# Transform the complete displacement vector EFCS to TDCS 
	(UxDn, UyDn, UzDn) =RotateObject3DNewCoords!(UxDn, UyDn, UzDn ,0.,0.,0.,VxR,VyR,VzR)
	(UxDss,UyDss,UzDss)=RotateObject3DNewCoords!(UxDss,UyDss,UzDss,0.,0.,0.,VxR,VyR,VzR)
	(UxDds,UyDds,UzDds)=RotateObject3DNewCoords!(UxDds,UyDds,UzDds,0.,0.,0.,VxR,VyR,VzR)

	#Flip these so we add component
	if P1[3]==0 && P2[3]==0 && P3[3]==0
		UzDn  = -UzDn;
		UzDss = -UzDss;
		UzDds = -UzDds;
	end

end

# Calculate first angular dislocation contribution POS
e13.=.-e13;
TDSetupD(x,y,z,A,Dn,Dss,Dds,nu,p1,e13, 
UxDn,UxDss,UxDds,
UyDn,UyDss,UyDds,
UzDn,UzDss,UzDds,Pos); 
e13.=.-e13; 
# Calculate second angular dislocation contribution
TDSetupD(x,y,z,B,Dn,Dss,Dds,nu,p2,e12,
UxDn,UxDss,UxDds,
UyDn,UyDss,UyDds,
UzDn,UzDss,UzDds,Pos); 
# Calculate third angular dislocation contribution
TDSetupD(x,y,z,C,Dn,Dss,Dds,nu,p3,e23,
UxDn,UxDss,UxDds,
UyDn,UyDss,UyDds,
UzDn,UzDss,UzDds,Pos);
 
# Calculate first angular dislocation contribution NEG
TDSetupD(x,y,z,A,Dn,Dss,Dds,nu,p1,e13,
UxDn,UxDss,UxDds,
UyDn,UyDss,UyDds,
UzDn,UzDss,UzDds,Neg);
# Calculate second angular dislocation contribution
e12.=.-e12;
TDSetupD(x,y,z,B,Dn,Dss,Dds,nu,p2,e12,
UxDn,UxDss,UxDds,
UyDn,UyDss,UyDds,
UzDn,UzDss,UzDds,Neg); e12.=.-e12;
# Calculate third angular dislocation contribution
e23.=.-e23;
TDSetupD(x,y,z,C,Dn,Dss,Dds,nu,p3,e23,
UxDn,UxDss,UxDds,
UyDn,UyDss,UyDds,
UzDn,UzDss,UzDds,Neg);	e23.=.-e23;
 
# Calculate the "incomplete" displacement vector components in TDCS
for i=eachindex(x)

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
(UxDn, UyDn, UzDn) =RotateObject3DNewCoords!(UxDn, UyDn, UzDn ,0.,0.,0.,Vnorm,Vstrike,Vdip)
(UxDss,UyDss,UzDss)=RotateObject3DNewCoords!(UxDss,UyDss,UzDss,0.,0.,0.,Vnorm,Vstrike,Vdip)
(UxDds,UyDds,UzDds)=RotateObject3DNewCoords!(UxDds,UyDds,UzDds,0.,0.,0.,Vnorm,Vstrike,Vdip)

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

end


function GlobalToTDECoords(P1,P2,P3,X,Y,Z,Vnorm,Vstrike,Vdip)
# Transform coordinates from EFCS into TDCS

#Inverse rot mat
Vx=zeros(3,1);Vy=zeros(3,1);Vz=zeros(3,1);
Vx[1]=Vnorm[1];  Vy[1]=Vnorm[2];  Vz[1]=Vnorm[3];
Vx[2]=Vstrike[1];Vy[2]=Vstrike[2];Vz[2]=Vstrike[3];
Vx[3]=Vdip[1];   Vy[3]=Vdip[2];   Vz[3]=Vdip[3];	


#Init some vars
p1 = copy(P1);
p2 = zeros(3,1);
p3 = copy(P3);
x=copy(X);y=copy(Y);z=copy(Z);

(x,y,z)=RotateObject3DNewCoords!(x,y,z,P2[1],P2[2],P2[3],Vx,Vy,Vz)
(p1[1],p1[2],p1[3])=RotateObject3DNewCoords!(p1[1],p1[2],p1[3],P2[1],P2[2],P2[3],Vx,Vy,Vz)
(p3[1],p3[2],p3[3])=RotateObject3DNewCoords!(p3[1],p3[2],p3[3],P2[1],P2[2],P2[3],Vx,Vy,Vz)#Vx,Vy,Vz

#Get interior angles and vectors along the triangle edges. 
(e12,e13,e23,A,B,C)=CalcTDVectsAndAngles(p1,p2,p3) #7 allocs
# Determine the best arteact-free configuration for each calculation point
(casepLog,casenLog,casezLog) = trimodefinder(y,z,x,p1[2:3],p2[2:3],p3[2:3]);
#Turn to cart index
casepLog=findall(casepLog);
casenLog=findall(casenLog)

return(p1,p2,p3,x,y,z,e12,e13,e23,A,B,C,casepLog,casenLog,casezLog)
end


function CalculateLocalTriCoords(P1,P2,P3)
# Calculate unit strike, dip and normal to TD vectors: For a horizontal TD 
# as an exception, if the normal vector points upward, the strike and dip 
# vectors point Northward and Westward, whereas if the normal vector points
# downward, the strike and dip vectors point Southward and Westward, 
# respectively.

eY = [0.;1.;0.];
eZ = [0.;0.;1.];

Vnorm = cross(P2-P1,P3-P1);
Vnorm = Vnorm/sqrt(Vnorm[1]^2+Vnorm[2]^2+Vnorm[3]^2)
Vstrike = cross(eZ,Vnorm);
if norm(Vstrike)==0
	Vstrike = eY*Vnorm[3];
end
Vstrike = Vstrike/sqrt(Vstrike[1]^2+Vstrike[2]^2+Vstrike[3]^2);
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
e12=e12/sqrt(e12[1]^2+e12[2]^2+e12[3]^2)  #could use :norm(e12);
e13=p3-p1;
e13=e13/sqrt(e13[1]^2+e13[2]^2+e13[3]^2);
e23=p3-p2;
e23=e23/sqrt(e23[1]^2+e23[2]^2+e23[3]^2);

# Calculate the TD angles
#A=e12'*e13;
A=(e12[1]*e13[1])+(e12[2]*e13[2])+(e12[3]*e13[3]);
A=acos(A);
B=(-e12[1]*e23[1])+(-e12[2]*e23[2])+(-e12[3]*e23[3]);
B=acos(B);
C=(e23[1]*e13[1])+(e23[2]*e13[2])+(e23[3]*e13[3]);
C=acos(C);

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

#Init some values outside of the loop
p2_2_m_p3_2=p2[2]-p3[2];
p3_1_m_p2_1=p3[1]-p2[1];
p1_1_m_p3_1=p1[1]-p3[1];
p1_2_m_p3_2=p1[2]-p3[2];
p3_2_m_p1_2=p3[2]-p1[2];
Base=(p2_2_m_p3_2*p1_1_m_p3_1+p3_1_m_p2_1*p1_2_m_p3_2);

trimode=ones(Int64, length(x),1)
for i=eachindex(x)

	a = (p2_2_m_p3_2*(x[i]-p3[1])+p3_1_m_p2_1*(y[i]-p3[2]))/Base;
	b = (p3_2_m_p1_2*(x[i]-p3[1])+p1_1_m_p3_1*(y[i]-p3[2]))/Base;
	c = 1.0-a-b;
		
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
for i=eachindex(trimode)
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
 UzDn,UzDss,UzDds,Index)
# TDSetupD transforms coordinates of the calculation points as well as 
# slip vector components from ADCS into TDCS. It then calculates the 
# displacements in ADCS and transforms them into TDCS.

Ct=SideVec[3];
St=SideVec[2];

#Here we compute slip vectors for a unit dislocation (of the given magnitude) (StrikeSlip and Dipslip separately)
# Transform the in-plane slip vector components from TDCS into ADCS
(Dss1,Dds0)=RotateObject2D(Dss,0.,0.,0.,Ct,St)
(Dss0,Dds1)=RotateObject2D(0.,Dds,0.,0.,Ct,St)


#Init arrays
Ang=-pi+alpha;
cosA = cos(Ang);
sinA = sin(Ang);

#Extra defs out of loop to speed it up
E1=(1.0-nu); #Elastic cons
E2=(1.0-2.0*nu);
cosA2=cosA^2;
sinADE1=sinA/8.0/pi/(1-nu);

Dn8p=Dn/8.0/pi;

 
# Calculate displacements associated with an angular dislocation in ADCS
for i=eachindex(Index) #1:length(x)
	
	Y=y[Index[i]];
	Z=z[Index[i]];
	
	# Transform coordinates of the calculation points from TDCS into ADCS
	(Y,Z)  =RotateObject2D!(Y,Z,TriVertex[2],TriVertex[3],Ct,St)
	
	(ux,uy,uz,vx,vy,vz,wx,wy,wz) = AngDisDisp(x[Index[i]],Y,Z,cosA,sinA,E1,E2,cosA2,sinADE1);

	
	if Dn!=0 #Only doing if needed
		#Ux is in global coords:
		#components due to opening
		UxDn[Index[i]]=UxDn[Index[i]]+(Dn8p/E1*ux);
		#Comp Uy and Uz in the current coords
		uyDn=Dn8p/E1*vx;
		uzDn=Dn8p/E1*wx;
		#Rotate to global coords			
		(uyDn,uzDn)=RotateObject2D!(uyDn,uzDn,0.,0.,Ct,-St)
		# #Add these to the total vector 
		UyDn[Index[i]]=UyDn[Index[i]]+uyDn;
		UzDn[Index[i]]=UzDn[Index[i]]+uzDn;		
	end
	
	#For Dss and Dds the local coordinates mean these must be combined 
	#to get the global contribution for these parts.
	
	if Dss!=0 #Only doing if needed
		#Ux is in global coords
		UxDss[Index[i]]=UxDss[Index[i]]+(Dss1/8.0/pi/E1*uy)+(Dds0*sinADE1*uz)
		#Comp Uy and Uz in the current coords 	
		uyDss=(Dss1*x[Index[i]]/8.0/pi/E1*vy)+(Dds0*x[Index[i]]*sinADE1*vz)	
		uzDss=(Dss1*x[Index[i]]/8.0/pi/E1*wy)+(Dds0*x[Index[i]]*sinADE1*wz)
		#Rotate to global coords			
		(uyDss,uzDss)=RotateObject2D!(uyDss,uzDss,0.,0.,Ct,-St) 
		#Add these to the total vector 		
		UyDss[Index[i]]=UyDss[Index[i]]+uyDss;
		UzDss[Index[i]]=UzDss[Index[i]]+uzDss;
	end		
		
	if Dds!=0 #Only doing if needed
		#Ux is in global coords
		UxDds[Index[i]]=UxDds[Index[i]]+(Dss0/8.0/pi/E1*uy)+(Dds1*sinADE1*uz)
		#Comp Uy and Uz in the current coords		
		uyDds=(Dss0*x[Index[i]]/8.0/pi/E1*vy)+(Dds1*x[Index[i]]*sinADE1*vz)		
		uzDds=(Dss0*x[Index[i]]/8.0/pi/E1*wy)+(Dds1*x[Index[i]]*sinADE1*wz)
		#Rotate to global coords		
		(uyDds,uzDds)=RotateObject2D!(uyDds,uzDds,0.,0.,Ct,-St) 
		#Add these to the total vector 		
		UyDds[Index[i]]=UyDds[Index[i]]+uyDds;
		UzDds[Index[i]]=UzDds[Index[i]]+uzDds;
	end

	
end

end



function AngDisDisp(x::Float64,y::Float64,z::Float64,
					cosA::Float64,sinA::Float64,
					E1::Float64,E2::Float64,
					cosA2::Float64,sinADE1::Float64)
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

#Commented parts on right are the orginal parts of Mehdi's function (burgers vector size). 
#This is now computed outside

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
return( ux::Float64,uy::Float64,uz::Float64,
		vx::Float64,vy::Float64,vz::Float64,
		wx::Float64,wy::Float64,wz::Float64)
end

function TDdisp_HarFunc(X,Y,Z,P1,P2,P3,Dss,Dds,Dn,nu,Vnorm,Vstrike,Vdip,
				UxDn, UyDn, UzDn,
				UxDss,UyDss,UzDss,
				UxDds,UyDds,UzDds)
# TDdisp_HarFunc calculates the harmonic function contribution to the
# displacements associated with a triangular dislocation in a half-space.
# The function cancels the surface normal tractions induced by the main and
# image dislocations.

# Calculate contribution of angular dislocation pair on each TD side 

# Side P1P2
AngSetupDispFSC(X,Y,Z,P1,P2,nu,Vnorm,Vstrike,Vdip,Dn,Dss,Dds,
UxDn, UyDn, UzDn,
UxDss,UyDss,UzDss,
UxDds,UyDds,UzDds); 
 
# Side P2P3
AngSetupDispFSC(X,Y,Z,P2,P3,nu,Vnorm,Vstrike,Vdip,Dn,Dss,Dds,
UxDn, UyDn, UzDn,
UxDss,UyDss,UzDss,
UxDds,UyDds,UzDds); 
 
# Side P3P1
AngSetupDispFSC(X,Y,Z,P3,P1,nu,Vnorm,Vstrike,Vdip,Dn,Dss,Dds,
UxDn, UyDn, UzDn,
UxDss,UyDss,UzDss,
UxDds,UyDds,UzDds); 

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

	#ALL ALLOCS BETWEEN HERE
	#println("reduce between here")
	###############################
    
	ey1=zeros(3);
	ey1[1:2] = SideVec[1:2];	
	ey1 = ey1/sqrt(ey1[1]^2+ey1[2]^2+ey1[3]^2)  #could use :norm(ey1);
	ey3 = -eZ;
	ey2 = cross(ey3,ey1);
	
	# Transform coordinates from EFCS to the second ADCS
	(y1B,y2B,y3B)=RotateObject3DNewCoords!(SideVec[1],SideVec[2],SideVec[3],0.,0.,0.,ey1,ey2,ey3)

	#Inverse rot mat
	VxR=zeros(3,1);VyR=zeros(3,1);VzR=zeros(3,1)
	VxR[1]=ey1[1];VyR[1]=ey1[2];VzR[1]=ey1[3];
	VxR[2]=ey2[1];VyR[2]=ey2[2];VzR[2]=ey2[3];
	VxR[3]=ey3[1];VyR[3]=ey3[2];VzR[3]=ey3[3];
	
	#println("Add the stuff above into loop, dont split the loop into two (Indx and indxf) instead just define beta in two ways at top.")
	###############################
	
	
	#Here we compute slip vectors for a unit dislocation (of the given magnitude) (Normal, StrikeSlip and Dipslip separately)	
	# Transform slip vector components from TDCS into EFCS
	(Dn1__,Dss0n_,Dds0n_) = RotateObject3DNewCoords(Dn,0.,0., 0.,0.,0.,Vnorm,Vstrike,Vdip);
	(Dn0ss,Dss1__,Dds1ss) = RotateObject3DNewCoords(0.,Dss,0.,0.,0.,0.,Vnorm,Vstrike,Vdip);
	(Dn0ds,Dss0ds,Dds1__) = RotateObject3DNewCoords(0.,0.,Dds,0.,0.,0.,Vnorm,Vstrike,Vdip);
	# Transform slip vector components from EFCS to ADCS
	(Dn1__,Dss0n_,Dds0n_)=RotateObject3DNewCoords!(Dn1__,Dss0n_,Dds0n_,0.,0.,0.,ey1,ey2,ey3)
	(Dn0ss,Dss1__,Dds1ss)=RotateObject3DNewCoords!(Dn0ss,Dss1__,Dds1ss,0.,0.,0.,ey1,ey2,ey3)
	(Dn0ds,Dss0ds,Dds1__)=RotateObject3DNewCoords!(Dn0ds,Dss0ds,Dds1__,0.,0.,0.,ey1,ey2,ey3);
	
	#indx=findall(I);
	pib=-pi+beta;
	sinPib = sin(pib);
	cosPib = cos(pib);
	cotPib = cot(pib);
	cotPib2=cot(pib/2);
	
	#Invert the bool
	#indxf=findall(.!I);
	b=beta;
	sinB = sin(b);
	cosB = cos(b);
	cotB = cot(b);
	cotB2=cot(b/2);
	
	# Configuration I
	for i=eachindex(X)
	
		x=X[i];
		y=Y[i];
		z=Z[i];
		# Transform coordinates from EFCS to the first ADCS
		(x,y,z)=RotateObject3DNewCoords!(x,y,z,PA[1],PA[2],PA[3],ey1,ey2,ey3)
		Bx=beta*x;
		if Bx>=0;
		#Call func that does the work
		#if I[i] == 1
			(uxA,uyA,uzA,vxA,vyA,vzA,wxA,wyA,wzA) = AngDisDispFSC(x,y,z,cosPib,sinPib,cotPib,cotPib2,nu,-PA[3]);
		else 
			(uxA,uyA,uzA,vxA,vyA,vzA,wxA,wyA,wzA) = AngDisDispFSC(x,y,z,cosB,sinB,cotB,cotB2,nu,-PA[3]);
		end
		
		#The local coordinates are such that the components must be combined 
		#to get the global contribution for these parts.
		if Dn!=0
			uxDn =-((Dn1__/4/pi/(1-nu)*uxA)+(Dss0n_/4/pi/(1-nu)*uyA)+(Dds0n_/4/pi/(1-nu)*uzA))
			uyDn =-((Dn1__/4/pi/(1-nu)*vxA)+(Dss0n_/4/pi/(1-nu)*vyA)+(Dds0n_/4/pi/(1-nu)*vzA))	
			uzDn =-((Dn1__/4/pi/(1-nu)*wxA)+(Dss0n_/4/pi/(1-nu)*wyA)+(Dds0n_/4/pi/(1-nu)*wzA))
		end			
		if Dss!=0		
			uxDss=-((Dn0ss/4/pi/(1-nu)*uxA)+(Dss1__/4/pi/(1-nu)*uyA)+(Dds1ss/4/pi/(1-nu)*uzA))
			uyDss=-((Dn0ss/4/pi/(1-nu)*vxA)+(Dss1__/4/pi/(1-nu)*vyA)+(Dds1ss/4/pi/(1-nu)*vzA))	
			uzDss=-((Dn0ss/4/pi/(1-nu)*wxA)+(Dss1__/4/pi/(1-nu)*wyA)+(Dds1ss/4/pi/(1-nu)*wzA))
		end
		if Dds!=0			
			uxDds=-((Dn0ds/4/pi/(1-nu)*uxA)+(Dss0ds/4/pi/(1-nu)*uyA)+(Dds1__/4/pi/(1-nu)*uzA))
			uyDds=-((Dn0ds/4/pi/(1-nu)*vxA)+(Dss0ds/4/pi/(1-nu)*vyA)+(Dds1__/4/pi/(1-nu)*vzA))
			uzDds=-((Dn0ds/4/pi/(1-nu)*wxA)+(Dss0ds/4/pi/(1-nu)*wyA)+(Dds1__/4/pi/(1-nu)*wzA))
		end			
		
		if Bx>=0;
		#Call func that does the work
		#if I[i] == 1		
			#Call func that does the work		
			(uxB,uyB,uzB,vxB,vyB,vzB,wxB,wyB,wzB) = AngDisDispFSC(x-y1B,y-y2B,z-y3B,cosPib,sinPib,cotPib,cotPib2,nu,-PB[3]);
		else
			(uxB,uyB,uzB,vxB,vyB,vzB,wxB,wyB,wzB) = AngDisDispFSC(x-y1B,y-y2B,z-y3B,cosB,sinB,cotB,cotB2,nu,-PB[3]);
		end

		if Dn!=0		
			#Add to the current value this part		
			uxDn  =uxDn  + ((Dn1__/4/pi/(1-nu)*uxB)+(Dss0n_/4/pi/(1-nu)*uyB)+(Dds0n_/4/pi/(1-nu)*uzB))
			uyDn  =uyDn  + ((Dn1__/4/pi/(1-nu)*vxB)+(Dss0n_/4/pi/(1-nu)*vyB)+(Dds0n_/4/pi/(1-nu)*vzB))		
			uzDn  =uzDn  + ((Dn1__/4/pi/(1-nu)*wxB)+(Dss0n_/4/pi/(1-nu)*wyB)+(Dds0n_/4/pi/(1-nu)*wzB))
			#Rotate to global coords	
			(uxDn, uyDn, uzDn) =RotateObject3DNewCoords!(uxDn, uyDn, uzDn ,0.,0.,0.,VxR,VyR,VzR)	
			#Add these to the total vector				
			UxDn[i]  =UxDn[i]  + uxDn;		
			UyDn[i]  =UyDn[i]  + uyDn;
			UzDn[i]  =UzDn[i]  + uzDn;			
		end
		if  Dss!=0
			#Add to the current value this part			
			uxDss =uxDss + ((Dn0ss/4/pi/(1-nu)*uxB)+(Dss1__/4/pi/(1-nu)*uyB)+(Dds1ss/4/pi/(1-nu)*uzB))
			uyDss =uyDss + ((Dn0ss/4/pi/(1-nu)*vxB)+(Dss1__/4/pi/(1-nu)*vyB)+(Dds1ss/4/pi/(1-nu)*vzB))	
			uzDss =uzDss + ((Dn0ss/4/pi/(1-nu)*wxB)+(Dss1__/4/pi/(1-nu)*wyB)+(Dds1ss/4/pi/(1-nu)*wzB))
			#Rotate to global coords		
			(uxDss,uyDss,uzDss)=RotateObject3DNewCoords!(uxDss,uyDss,uzDss,0.,0.,0.,VxR,VyR,VzR)	
			#Add these to the total vector				
			UxDss[i] =UxDss[i] + uxDss;
			UyDss[i] =UyDss[i] + uyDss;		
			UzDss[i] =UzDss[i] + uzDss;			
		end
		if  Dds!=0		
			#Add to the current value this part			
			uxDds =uxDds + ((Dn0ds/4/pi/(1-nu)*uxB)+(Dss0ds/4/pi/(1-nu)*uyB)+(Dds1__/4/pi/(1-nu)*uzB))
			uyDds =uyDds + ((Dn0ds/4/pi/(1-nu)*vxB)+(Dss0ds/4/pi/(1-nu)*vyB)+(Dds1__/4/pi/(1-nu)*vzB))
			uzDds =uzDds + ((Dn0ds/4/pi/(1-nu)*wxB)+(Dss0ds/4/pi/(1-nu)*wyB)+(Dds1__/4/pi/(1-nu)*wzB))
			#Rotate to global coords			
			(uxDds,uyDds,uzDds)=RotateObject3DNewCoords!(uxDds,uyDds,uzDds,0.,0.,0.,VxR,VyR,VzR)
			#Add these to the total vector			
			UxDds[i] =UxDds[i] + uxDds;
			UyDds[i] =UyDds[i] + uyDds;
			UzDds[i] =UzDds[i] + uzDds;			
		end			
		
	end

end	#if statement

end


function CalcSideVec(PA,PB)
# Calculate TD side vector and the angle of the angular dislocation pair

SideVec = PB-PA;
eZ = [0.;0.;1.];


G=-SideVec'*eZ/sqrt(SideVec[1]^2+SideVec[2]^2+SideVec[3]^2) 
beta = acos(G[1]);

return(SideVec,eZ,beta)
end


function AngDisDispFSC( y1::Float64,y2::Float64,y3::Float64,
						cosB::Float64,sinB::Float64,cotB::Float64,cotB2::Float64,
						nu::Float64,a::Float64)
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
return( ux::Float64,uy::Float64,uz::Float64,
		vx::Float64,vy::Float64,vz::Float64,
		wx::Float64,wy::Float64,wz::Float64)

end


function TDstrainFS(X,Y,Z,P1,P2,P3,Dss,Dds,Dn,mu,lambda,nu,ImageFlag,Vnorm,Vstrike,Vdip,
					p1,p2,p3,x,y,z,e12,e13,e23,A,B,C,Pos,Neg,casezLog,
		 ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
		 ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
		 ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds)
# TDstressFS 
# Calculates stresses and strains associated with a triangular dislocation 
# in an elastic full-space.


if ImageFlag==1; #This means we are computing the iamge dislocation
	RotInvMat=zeros(1,9); #Single alloc here!
	RotInvMat[1]=Vnorm[1]; RotInvMat[2]=Vstrike[1]; RotInvMat[3]=Vdip[1];
	RotInvMat[4]=Vnorm[2]; RotInvMat[5]=Vstrike[2]; RotInvMat[6]=Vdip[2];
	RotInvMat[7]=Vnorm[3]; RotInvMat[8]=Vstrike[3]; RotInvMat[9]=Vdip[3];
	# Transform the strain tensor components from TDCS into EFCS
	(ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn) = TensorTransformation3D!(ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn,RotInvMat);
	(ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss) = TensorTransformation3D!(ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,RotInvMat);
	(ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds) = TensorTransformation3D!(ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,RotInvMat);
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

# Calculate first angular dislocation contribution POS
e13.=.-e13;
TDSetupS(x,y,z,A,Dn,Dss,Dds,nu,p1,e13,
ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,Pos); e13.=.-e13; #3allocs
 
# Calculate second angular dislocation contribution
TDSetupS(x,y,z,B,Dn,Dss,Dds,nu,p2,e12,
ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,Pos); 
 
# Calculate third angular dislocation contribution
TDSetupS(x,y,z,C,Dn,Dss,Dds,nu,p3,e23,
ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,Pos);


# Calculate first angular dislocation contribution NEG
TDSetupS(x,y,z,A,Dn,Dss,Dds,nu,p1,e13,
ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,Neg);

# Calculate second angular dislocation contribution
e12.=.-e12;
TDSetupS(x,y,z,B,Dn,Dss,Dds,nu,p2,e12,
ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,Neg); e12.=.-e12;
 
# Calculate third angular dislocation contribution 
e23.=.-e23;
TDSetupS(x,y,z,C,Dn,Dss,Dds,nu,p3,e23,
ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,Neg); e23.=.-e23;
		
# Calculate the strain tensor components in TDCS
for i=eachindex(x)
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

RotMat=zeros(1,9); #Single alloc here!
RotMat[1:3]=Vnorm; RotMat[4:6]=Vstrike; RotMat[7:9]=Vdip; 
# Transform the strain tensor components from TDCS into EFCS
(ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn) = TensorTransformation3D!(ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn,RotMat);
(ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss) = TensorTransformation3D!(ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,RotMat);
(ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds) = TensorTransformation3D!(ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,RotMat);

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
 
end


function TDSetupS(x,y,z,alpha,Dn,Dss,Dds,nu,TriVertex,SideVec,
				 ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
				 ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
				 ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,Index)
# TDSetupS transforms coordinates of the calculation points as well as 
# slip vector components from ADCS into TDCS. It then calculates the 
# strains in ADCS and transforms them into TDCS.

Ct=SideVec[3];
St=SideVec[2];

#Here we compute slip vectors for a unit dislocation (of the given magnitude) (StrikeSlip and Dipslip separately)	
# Transform the in-plane slip vector components from TDCS into ADCS
(Dss1,Dds0)=RotateObject2D(Dss,0.,0.,0.,Ct,St)
(Dss0,Dds1)=RotateObject2D(0.,Dds,0.,0.,Ct,St)

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
#B=[[1. 0. 0.];[0. Ct St];[0. -St Ct]]	 # 3x3 Transformation matrix
B=zeros(1,9); #Single alloc here!
B[1]=1.; B[5]=Ct; B[6]=-St; B[8]=St; B[9]=Ct; #Col wise indexing

# Calculate strains associated with an angular dislocation in ADCS
for i=eachindex(Index)

	Y=y[Index[i]];
	Z=z[Index[i]];
	
	# Transform coordinates of the calculation points from TDCS into ADCS
	(Y,Z)  =RotateObject2D!(Y,Z,TriVertex[2],TriVertex[3],Ct,St)

	(exxDn,exxDss,exxDds,
	 eyyDn,eyyDss,eyyDds,
	 ezzDn,ezzDss,ezzDds,
	 exyDn,exyDss,exyDds,
	 exzDn,exzDss,exzDds,
	 eyzDn,eyzDss,eyzDds,
	 rFi_rx,rFi_ry,rFi_rz) =
	 AngDisStrain(x[Index[i]],Y,Z,cosA,sinA,nu,E2,sinADE1)

	if Dn!=0 #Only doing if needed
		exxdn= Dn*(rFi_rx)+E3* exxDn
		eyydn= E3* eyyDn; 
		ezzdn= E3* ezzDn;
		exydn= Dn*(rFi_ry)/2-E3* exyDn; 
		exzdn= Dn*(rFi_rz)/2-E3* exzDn;
		eyzdn= E3* eyzDn;
		#Rotate to global coords
		(exxdn,eyydn,ezzdn,exydn,exzdn,eyzdn) = 
		TensorTransformation3D!(exxdn,eyydn,ezzdn,exydn,exzdn,eyzdn,B);
		#Add to total vector
		ExxDn[Index[i]] 	= ExxDn[Index[i]] +exxdn;
		EyyDn[Index[i]] 	= EyyDn[Index[i]] +eyydn
		EzzDn[Index[i]] 	= EzzDn[Index[i]] +ezzdn
		ExyDn[Index[i]] 	= ExyDn[Index[i]] +exydn
		ExzDn[Index[i]] 	= ExzDn[Index[i]] +exzdn
		EyzDn[Index[i]] 	= EyzDn[Index[i]] +eyzdn
	end		
	
	#For Dss and Dds the local coordinates mean these must be combined 
	#to get the global contribution for these parts.
	if Dss!=0 #Only doing if needed
		#A constant
		E4Dss1=Dss1*x[Index[i]]/8/pi/E1;	
		#Comp mixed components (in the current coords)
		exxdss= (-E4Dss1* exxDss)  +  (Dds0*x[Index[i]]*sinADE1* exxDds);
		eyydss= (Dss1*(rFi_ry)-E4Dss1* eyyDss)  +  (Dds0*x[Index[i]]*sinADE1* eyyDds);
		ezzdss= (-E4Dss1* ezzDss) + (Dds0*(rFi_rz)+Dds0*x[Index[i]]*sinADE1* ezzDds); 
		exydss= (Dss1*(rFi_rx)/2+E5Dss1* exyDss)  -  (Dds0*sinADE1* exyDds);
		exzdss= (E5Dss1* exzDss) + (Dds0*(rFi_rx)/2-Dds0*sinADE1*exzDds) ; 
		eyzdss= (Dss1*(rFi_rz)/2-E4Dss1* eyzDss) + (Dds0*(rFi_ry)/2-Dds0*x[Index[i]]*sinADE1* eyzDds);
		#Rotate to global coords		
		(exxdss,eyydss,ezzdss,exydss,exzdss,eyzdss) = 
		TensorTransformation3D!(exxdss,eyydss,ezzdss,exydss,exzdss,eyzdss,B);	
		#Add to total vector		
		ExxDss[Index[i]] 	= ExxDss[Index[i]] +exxdss;
		EyyDss[Index[i]] 	= EyyDss[Index[i]] +eyydss
		EzzDss[Index[i]] 	= EzzDss[Index[i]] +ezzdss
		ExyDss[Index[i]] 	= ExyDss[Index[i]] +exydss
		ExzDss[Index[i]] 	= ExzDss[Index[i]] +exzdss
		EyzDss[Index[i]] 	= EyzDss[Index[i]] +eyzdss
	end		
		
	if Dds!=0 #Only doing if needed	
		#A constant
		E4Dss0=Dss0*x[Index[i]]/8/pi/E1;
		#Comp mixed components (in the current coords)		
		exxdds= (-E4Dss0* exxDss)  +  (Dds1*x[Index[i]]*sinADE1* exxDds);
		eyydds= (Dss0*(rFi_ry)-E4Dss0* eyyDss)  +  (Dds1*x[Index[i]]*sinADE1* eyyDds);
		ezzdds= (-E4Dss0* ezzDss) + (Dds1*(rFi_rz)+Dds1*x[Index[i]]*sinADE1* ezzDds); 
		exydds= (Dss0*(rFi_rx)/2+E5Dss0* exyDss)  -  (Dds1*sinADE1* exyDds);	
		exzdds= (E5Dss0* exzDss) + (Dds1*(rFi_rx)/2-Dds1*sinADE1*exzDds);
		eyzdds= (Dss0*(rFi_rz)/2-E4Dss0* eyzDss) + (Dds1*(rFi_ry)/2-Dds1*x[Index[i]]*sinADE1* eyzDds);		  
		#Rotate to global coords
		(exxdds,eyydds,ezzdds,exydds,exzdds,eyzdds) = 
		TensorTransformation3D!(exxdds,eyydds,ezzdds,exydds,exzdds,eyzdds,B);
		#Add to total vector		
		ExxDds[Index[i]] 	= ExxDds[Index[i]] +exxdds
		EyyDds[Index[i]] 	= EyyDds[Index[i]] +eyydds
		EzzDds[Index[i]] 	= EzzDds[Index[i]] +ezzdds
		ExyDds[Index[i]] 	= ExyDds[Index[i]] +exydds
		ExzDds[Index[i]] 	= ExzDds[Index[i]] +exzdds
		EyzDds[Index[i]] 	= EyzDds[Index[i]] +eyzdds
	end		
end	

end


function AngDisStrain(x::Float64,y::Float64,z::Float64,
					  cosA::Float64,sinA::Float64,
					  nu::Float64,E2::Float64,sinADE1::Float64)
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

return(ExxBx::Float64,ExxBy::Float64,ExxBz::Float64,
	   EyyBx::Float64,EyyBy::Float64,EyyBz::Float64,
	   EzzBx::Float64,EzzBy::Float64,EzzBz::Float64,
	   ExyBx::Float64,ExyBy::Float64,ExyBz::Float64,
	   ExzBx::Float64,ExzBy::Float64,ExzBz::Float64,
	   EyzBx::Float64,EyzBy::Float64,EyzBz::Float64,
	   rFi_rx::Float64,rFi_ry::Float64,rFi_rz::Float64)
end


function TDstrain_HarFunc(X,Y,Z,P1,P2,P3,Dss,Dds,Dn,mu,lambda,nu,Vnorm,Vstrike,Vdip,
			 ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
			 ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
			 ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds);
# TDstrain_HarFunc calculates the harmonic function contribution to the
# strains and stresses associated with a triangular dislocation in a 
# half-space. The function cancels the surface normal tractions induced by 
# the main and image dislocations.

# Calculate contribution of angular dislocation pair on each TD side 
# P1P2
AngSetupStrainFSC(X,Y,Z,Dn,Dss,Dds,P1,P2,mu,lambda,nu,Vnorm,Vstrike,Vdip,
ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds);
 
# P2P3
AngSetupStrainFSC(X,Y,Z,Dn,Dss,Dds,P2,P3,mu,lambda,nu,Vnorm,Vstrike,Vdip,
ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds); 
 
 
# P3P1 
AngSetupStrainFSC(X,Y,Z,Dn,Dss,Dds,P3,P1,mu,lambda,nu,Vnorm,Vstrike,Vdip,
ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
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
    
	ey1=zeros(3);
	ey1[1:2] = SideVec[1:2];	
	ey1 = ey1/sqrt(ey1[1]^2+ey1[2]^2+ey1[3]^2)  #could use :norm(ey1);
	ey3 = -eZ;
	ey2 = cross(ey3,ey1);
	


	# Transform coordinates from EFCS to the second ADCS
	(y1B,y2B,y3B)=RotateObject3DNewCoords!(SideVec[1],SideVec[2],SideVec[3],0.,0.,0.,ey1,ey2,ey3)

	#Inverse rot mat
	AFlip=zeros(1,9);
	AFlip[1]=ey1[1];AFlip[2]=ey1[2];AFlip[3]=ey1[3];
	AFlip[4]=ey2[1];AFlip[5]=ey2[2];AFlip[6]=ey2[3];
	AFlip[7]=ey3[1];AFlip[8]=ey3[2];AFlip[9]=ey3[3];
	#AFlip = [ey1[1] ey1[2] ey1[3] ey2[1] ey2[2] ey2[3] ey3[1] ey3[2] ey3[3]]; # Transformation matri
	
	#Inverse rot mat
	VxR=zeros(3,1);VyR=zeros(3,1);VzR=zeros(3,1)
	VxR[1]=ey1[1];VyR[1]=ey1[2];VzR[1]=ey1[3];
	VxR[2]=ey2[1];VyR[2]=ey2[2];VzR[2]=ey2[3];
	VxR[3]=ey3[1];VyR[3]=ey3[2];VzR[3]=ey3[3];
	
	
	#Here we compute slip vectors for a unit dislocation (of the given magnitude) (Normal, StrikeSlip and Dipslip separately)		
	## Transform slip vector components from TDCS into EFCS
	(Dn1__i,Dss0n_i,Dds0n_) = RotateObject3DNewCoords(Dn,0.,0., 0.,0.,0.,Vnorm,Vstrike,Vdip);
	(Dn0ssi,Dss1__i,Dds1ss) = RotateObject3DNewCoords(0.,Dss,0.,0.,0.,0.,Vnorm,Vstrike,Vdip);
	(Dn0dsi,Dss0dsi,Dds1__) = RotateObject3DNewCoords(0.,0.,Dds,0.,0.,0.,Vnorm,Vstrike,Vdip);
	# Transform slip vector components from EFCS to ADCS
	(Dn1__i,Dss0n_i,Dds0n_)=RotateObject3DNewCoords(Dn1__i,Dss0n_i,Dds0n_,0.,0.,0.,ey1,ey2,ey3)
	(Dn0ssi,Dss1__i,Dds1ss)=RotateObject3DNewCoords(Dn0ssi,Dss1__i,Dds1ss,0.,0.,0.,ey1,ey2,ey3)
	(Dn0dsi,Dss0dsi,Dds1__)=RotateObject3DNewCoords(Dn0dsi,Dss0dsi,Dds1__,0.,0.,0.,ey1,ey2,ey3)
    
	#Init some vars
	#indx=findall(I);
	pib=pi-beta;
	sinPib = sin(pib);
	cosPib = cos(pib);
	cotPib = cot(pib);
	
	#Invert the bool
	#indxf=findall(.!I);
	b=beta;
	sinB = sin(b);
	cosB = cos(b);
	cotB = cot(b);
	
    # Configuration I
	for i=eachindex(X)
	
		x=X[i];
		y=Y[i];
		z=Z[i];
		# Transform coordinates from EFCS to the first ADCS
		(x,y,z)=RotateObject3DNewCoords!(x,y,z,PA[1],PA[2],PA[3],ey1,ey2,ey3)
		Bx=beta*x;
		if Bx>=0;
		#Call func that does the work		
		#if I[i] == 1		
			 #xx    yy    zz    xy    xz    yz
			(v11Ax,v22Ax,v33Ax,v12Ax,v13Ax,v23Ax,  #Dn
			 v11Ay,v22Ay,v33Ay,v12Ay,v13Ay,v23Ay,  #Dss
			 v11Az,v22Az,v33Az,v12Az,v13Az,v23Az) = AngDisStrainFSC(-x,-y,z,cosPib,sinPib,cotPib,nu,-PA[3]);
			v13Ax = -v13Ax;
			v13Ay = -v13Ay;
			v13Az = -v13Az;
			v23Ax = -v23Ax;
			v23Ay = -v23Ay;
			v23Az = -v23Az;
			Dn1__=-Dn1__i;
			Dn0ss=-Dn0ssi;
			Dn0ds=-Dn0dsi;
			Dss0n_=-Dss0n_i;
			Dss1__=-Dss1__i;
			Dss0ds=-Dss0dsi;
		else
			 #xx    yy    zz    xy    xz    yz
			(v11Ax,v22Ax,v33Ax,v12Ax,v13Ax,v23Ax,  #Dn
			 v11Ay,v22Ay,v33Ay,v12Ay,v13Ay,v23Ay,  #Dss
			 v11Az,v22Az,v33Az,v12Az,v13Az,v23Az) = AngDisStrainFSC(x,y,z,cosB,sinB,cotB,nu,-PA[3]);		
			Dn1__=Dn1__i;
			Dn0ss=Dn0ssi;
			Dn0ds=Dn0dsi;
			Dss0n_=Dss0n_i;
			Dss1__=Dss1__i;
			Dss0ds=Dss0dsi;		
		end
		 
		#The local coordinates are such that the components must be combined 
		#to get the global contribution for these parts.
		if Dn!=0	
			exxDn = -((Dn1__*v11Ax)+(Dss0n_*v11Ay)+(Dds0n_*v11Az))
			eyyDn = -((Dn1__*v22Ax)+(Dss0n_*v22Ay)+(Dds0n_*v22Az))
			ezzDn = -((Dn1__*v33Ax)+(Dss0n_*v33Ay)+(Dds0n_*v33Az))
			exyDn = -((Dn1__/2*v12Ax)+(Dss0n_/2*v12Ay)+(Dds0n_/2*v12Az))
			exzDn = -((Dn1__/2*v13Ax)+(Dss0n_/2*v13Ay)+(Dds0n_/2*v13Az))
			eyzDn = -((Dn1__/2*v23Ax)+(Dss0n_/2*v23Ay)+(Dds0n_/2*v23Az))
		end
		if Dss!=0			
			exxDss= -((Dn0ss*v11Ax)+(Dss1__*v11Ay)+(Dds1ss*v11Az))
			eyyDss= -((Dn0ss*v22Ax)+(Dss1__*v22Ay)+(Dds1ss*v22Az))
			ezzDss= -((Dn0ss*v33Ax)+(Dss1__*v33Ay)+(Dds1ss*v33Az))
			exyDss= -((Dn0ss/2*v12Ax)+(Dss1__/2*v12Ay)+(Dds1ss/2*v12Az))
			exzDss= -((Dn0ss/2*v13Ax)+(Dss1__/2*v13Ay)+(Dds1ss/2*v13Az))
			eyzDss= -((Dn0ss/2*v23Ax)+(Dss1__/2*v23Ay)+(Dds1ss/2*v23Az))
		end
		if Dds!=0		
			exxDds= -((Dn0ds*v11Ax)+(Dss0ds*v11Ay)+(Dds1__*v11Az))		
			eyyDds= -((Dn0ds*v22Ax)+(Dss0ds*v22Ay)+(Dds1__*v22Az))
			ezzDds= -((Dn0ds*v33Ax)+(Dss0ds*v33Ay)+(Dds1__*v33Az))
			exyDds= -((Dn0ds/2*v12Ax)+(Dss0ds/2*v12Ay)+(Dds1__/2*v12Az))
			exzDds= -((Dn0ds/2*v13Ax)+(Dss0ds/2*v13Ay)+(Dds1__/2*v13Az))				
			eyzDds= -((Dn0ds/2*v23Ax)+(Dss0ds/2*v23Ay)+(Dds1__/2*v23Az))
		end		

		if Bx>=0;
		#Call func that does the work		
		#if I[i] == 1
			(v11Bx,v22Bx,v33Bx,v12Bx,v13Bx,v23Bx,
			 v11By,v22By,v33By,v12By,v13By,v23By,
			 v11Bz,v22Bz,v33Bz,v12Bz,v13Bz,v23Bz) = AngDisStrainFSC(-x+y1B,-y+y2B,z-y3B,cosPib,sinPib,cotPib,nu,-PB[3]);
			v13Bx = -v13Bx;
			v13By = -v13By;
			v13Bz = -v13Bz;
			v23Bx = -v23Bx;
			v23By = -v23By;
			v23Bz = -v23Bz;			
		else
			(v11Bx,v22Bx,v33Bx,v12Bx,v13Bx,v23Bx,
			 v11By,v22By,v33By,v12By,v13By,v23By,
			 v11Bz,v22Bz,v33Bz,v12Bz,v13Bz,v23Bz) = AngDisStrainFSC(x-y1B,y-y2B,z-y3B,cosB,sinB,cotB,nu,-PB[3]);
		end			 
		 
		if Dn!=0		
			#Add to the current value this part			
			exxDn = exxDn + ((Dn1__*v11Bx)+(Dss0n_*v11By)+(Dds0n_*v11Bz))
			eyyDn = eyyDn + ((Dn1__*v22Bx)+(Dss0n_*v22By)+(Dds0n_*v22Bz))
			ezzDn = ezzDn + ((Dn1__*v33Bx)+(Dss0n_*v33By)+(Dds0n_*v33Bz))
			exyDn = exyDn + ((Dn1__/2*v12Bx)+(Dss0n_/2*v12By)+(Dds0n_/2*v12Bz))
			exzDn = exzDn + ((Dn1__/2*v13Bx)+(Dss0n_/2*v13By)+(Dds0n_/2*v13Bz))
			eyzDn = eyzDn + ((Dn1__/2*v23Bx)+(Dss0n_/2*v23By)+(Dds0n_/2*v23Bz))
			#transform to global coords				
			(exxDn,eyyDn,ezzDn,exyDn,exzDn,eyzDn) = TensorTransformation3D!(exxDn,eyyDn,ezzDn,exyDn,exzDn,eyzDn,AFlip);
			#Add these to the total vector		
			ExxDn[i] = ExxDn[i] +exxDn;
			EyyDn[i] = EyyDn[i] +eyyDn;
			EzzDn[i] = EzzDn[i] +ezzDn;
			ExyDn[i] = ExyDn[i] +exyDn;
			ExzDn[i] = ExzDn[i] +exzDn;
			EyzDn[i] = EyzDn[i] +eyzDn;		
		end			
		if Dss!=0		
			#Add to the current value this part			
			exxDss= exxDss+ ((Dn0ss*v11Bx)+(Dss1__*v11By)+(Dds1ss*v11Bz))
			eyyDss= eyyDss+ ((Dn0ss*v22Bx)+(Dss1__*v22By)+(Dds1ss*v22Bz))
			ezzDss= ezzDss+ ((Dn0ss*v33Bx)+(Dss1__*v33By)+(Dds1ss*v33Bz))
			exyDss= exyDss+ ((Dn0ss/2*v12Bx)+(Dss1__/2*v12By)+(Dds1ss/2*v12Bz))
			exzDss= exzDss+ ((Dn0ss/2*v13Bx)+(Dss1__/2*v13By)+(Dds1ss/2*v13Bz))
			eyzDss= eyzDss+ ((Dn0ss/2*v23Bx)+(Dss1__/2*v23By)+(Dds1ss/2*v23Bz))
			#transform to global coords	
			(exxDss,eyyDss,ezzDss,exyDss,exzDss,eyzDss) = TensorTransformation3D!(exxDss,eyyDss,ezzDss,exyDss,exzDss,eyzDss,AFlip);
			#Add these to the total vector	
			ExxDss[i]= ExxDss[i]+exxDss;
			EyyDss[i]= EyyDss[i]+eyyDss;	
			EzzDss[i]= EzzDss[i]+ezzDss;
			ExyDss[i]= ExyDss[i]+exyDss;		
			ExzDss[i]= ExzDss[i]+exzDss;	
			EyzDss[i]= EyzDss[i]+eyzDss;			
		end
		if Dds!=0	
			#Add to the current value this part			
			exxDds= exxDds+ ((Dn0ds*v11Bx)+(Dss0ds*v11By)+(Dds1__*v11Bz))		
			eyyDds= eyyDds+ ((Dn0ds*v22Bx)+(Dss0ds*v22By)+(Dds1__*v22Bz))		
			ezzDds= ezzDds+ ((Dn0ds*v33Bx)+(Dss0ds*v33By)+(Dds1__*v33Bz))	
			exyDds= exyDds+ ((Dn0ds/2*v12Bx)+(Dss0ds/2*v12By)+(Dds1__/2*v12Bz))	
			exzDds= exzDds+ ((Dn0ds/2*v13Bx)+(Dss0ds/2*v13By)+(Dds1__/2*v13Bz))	
			eyzDds= eyzDds+ ((Dn0ds/2*v23Bx)+(Dss0ds/2*v23By)+(Dds1__/2*v23Bz))	
			#transform to global coords			
			(exxDds,eyyDds,ezzDds,exyDds,exzDds,eyzDds) = TensorTransformation3D!(exxDds,eyyDds,ezzDds,exyDds,exzDds,eyzDds,AFlip);	
			#Add these to the total vector				
			ExxDds[i]= ExxDds[i]+exxDds;	
			EyyDds[i]= EyyDds[i]+eyyDds;
			EzzDds[i]= EzzDds[i]+ezzDds;
			ExyDds[i]= ExyDds[i]+exyDds;	
			ExzDds[i]= ExzDds[i]+exzDds;		
			EyzDds[i]= EyzDds[i]+eyzDds;	
		end
	end
    
		
end


end


function AngDisStrainFSC(y1::Float64,y2::Float64,y3::Float64,
						 cosB::Float64,sinB::Float64,cotB::Float64,
						 nu::Float64,a::Float64)
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
	

return(v11x::Float64,v22x::Float64,v33x::Float64,v12x::Float64,v13x::Float64,v23x::Float64,
	   v11y::Float64,v22y::Float64,v33y::Float64,v12y::Float64,v13y::Float64,v23y::Float64,
	   v11z::Float64,v22z::Float64,v33z::Float64,v12z::Float64,v13z::Float64,v23z::Float64)
end