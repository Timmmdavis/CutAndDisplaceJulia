function TD(X,Y,Z,P1,P2,P3,Dss,Dds,Dn,nu,mu,DispFlag,StressFlag,HSflag)
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


bx = Dn;  # DisplacementNormal
by = Dss; # DisplacementStrike-slip
bz = Dds; # DisplacementDip-slip


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
	(ue,un,uv) = TDdispFS(X,Y,Z,P1,P2,P3,by,bz,bx,nu);
	
	if HSflag==1
		
		# Calculate harmonic function contribution to displacements
		(ueFSC,unFSC,uvFSC) = TDdisp_HarFunc(X,Y,Z,P1,P2,P3,by,bz,bx,nu);

		# Calculate image dislocation contribution to displacements
		(ueIS,unIS,uvIS) = TDdispFS(X,Y,Z,P1i,P2i,P3i,by,bz,bx,nu);
		if P1i[3]==0 && P2i[3]==0 && P3i[3]==0
			uvIS = -uvIS;
		end

		# Calculate the complete displacement vector components in EFCS
		ue = ue+ueIS+ueFSC;
		un = un+unIS+unFSC;
		uv = uv+uvIS+uvFSC;

		if P1i[3]==0 && P2i[3]==0 && P3i[3]==0
			ue = -ue;
			un = -un;
			uv = -uv;
		end
		
	end
	
end 

if StressFlag==1

	#Elastic con
	lambda=(2*mu*nu)/(1-(2*nu));

	(Exx,Eyy,Ezz,Exy,Exz,Eyz)=TDstrainFS(X,Y,Z,P1,P2,P3,by,bz,bx,mu,lambda)	
	
	if HSflag==1
	
		(ExxFSC,EyyFSC,EzzFSC,ExyFSC,ExzFSC,EyzFSC) = TDstrain_HarFunc(X,Y,Z,P1,P2,P3,by,bz,bx,mu,lambda,nu);
	
		# Calculate image dislocation contribution to strains and stresses
		(ExxIS,EyyIS,EzzIS,ExyIS,ExzIS,EyzIS) = TDstrainFS(X,Y,Z,P1i,P2i,P3i,by,bz,bx,mu,lambda);
		if P1i[3]==0 && P2i[3]==0 && P3i[3]==0
			ExzIS = -ExzIS;
			EyzIS = -EyzIS;
		end
		
		Exx=Exx+ExxFSC+ExxIS;
		Eyy=Eyy+EyyFSC+EyyIS;
		Ezz=Ezz+EzzFSC+EzzIS;
		Exy=Exy+ExyFSC+ExyIS;
		Exz=Exz+ExzFSC+ExzIS;
		Eyz=Eyz+EyzFSC+EyzIS;
		
	end
	
end

if DispFlag==0
	return(Exx,Eyy,Ezz,Exy,Exz,Eyz)
	
elseif StressFlag==0
	return(ue,un,uv)
	
else
	return(Exx,Eyy,Ezz,Exy,Exz,Eyz,ue,un,uv) #sum outside when needed
	
end

end






