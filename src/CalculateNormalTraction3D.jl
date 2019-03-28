function CalculateNormalTraction3D( Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,CosAx,CosAy,CosAz )
# CalculateNormalTraction3d: Calculates normal traction on a given 3D plane
#                   with Cartesian stress tensor components defined.
#                   Based on Equation 6.49 in Pollard, D.D. and Fletcher,
#                   R.C., 2005. Fundamentals of structural geology.
#                   Cambridge University Press.
#               
# usage #1:
# [ NormalTraction ] = CalculateNormalTraction3d( Pxx,Pyy,Pzz,Pxy,Pxz,Pyz,CosAx,CosAy,CosAz )
#
# Arguments: (input)
# Pxx,Pyy,Pzz 
# Pxy,Pxz,Pyz 		- The stress tensor components on this plane.
#                    (column vectors)
#
# CosAx,CosAy,CosAz - The seperated direction cosines of the surface 
#                    (column vectors)
#
# Arguments: (output)
# NormalTraction    - Column vector of normal traction magnitudes on this
#                     plane. Notation in scripts is typically: Tn.
#
# Example usage 1:
#
# # Get normal traction for a plane dipping 45 degrees facing east:
# CosAx=cosd(45);
# CosAy=0;
# CosAz=cosd(45);
# Sxx=1; Syy=0; Szz=0; Sxy=0; Sxz=0; Syz=0;
# [ NormalTraction ] = CalculateNormalTraction3d( Sxx,Syy,Szz,Sxy,Sxz,Syz,CosAx,CosAy,CosAz )
#
#  Author: Tim Davis
#  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

#Init vector to fill
NormalTraction=zeros(size(Pxx));
for i=1:size(Pxx,1) #On each col

	CosAx2=CosAx[i].^2;
	CosAy2=CosAy[i].^2;
	CosAz2=CosAz[i].^2;
	CosAxy=CosAx[i]*CosAy[i];
	CosAyz=CosAy[i]*CosAz[i];
	CosAzx=CosAz[i]*CosAx[i];
	
	for j=1:size(Pxx,2) #Running through the rows..
		# Normal traction on the planes. 
		# Equation 6.49 Pollard and Fletcher Fundamentals
		AA=Pxx[i,j]*CosAx2;
		BB=Pyy[i,j]*CosAy2;
		CC=Pzz[i,j]*CosAz2;
		DD=2*Pxy[i,j]*CosAxy;
		EE=2*Pyz[i,j]*CosAyz;
		FF=2*Pxz[i,j]*CosAzx;
		NormalTraction[i,j]=AA+BB+CC+DD+EE+FF;
	end
	
end
return NormalTraction
end
