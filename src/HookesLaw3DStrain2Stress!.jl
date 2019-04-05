function HookesLaw3DStrain2Stress!( Exx,Eyy,Ezz,Exy,Exz,Eyz,λ,G )
# HookesLaw3dStrain2Stress: Using 3D Hooke's law
#                   Converting the strain tensors to stress tensors using
#                   the Young's modulus, shear modulus and Poisson's ratio.
#                   Works with column vectors or single values.
#
#                   Eq. 7.131 and 7.132 of Pollard and Fletcher, 2005, 
#
#               
# usage #1:
# [ Sxx,Syy,Szz,Sxy,Sxz,Syz ] = HookesLaw3dStrain2Stress( Exx,Eyy,Ezz,Exy,Exz,Eyz,λ,G )
#
# Arguments: (input)
# Exx,Eyy,Ezz       
# Exy,Exz,Eyz       - Strain tensor components. (Can be column vectors). 
#
# λ            - Lame's constant
#
# G                - The shear modulus
#
# Arguments: (output)
# Sxx,Syy,Szz       
# Sxy,Sxz,Syz        - Stress tensor components. 
#
# Example usage:
#
# G=5; #[Gpa]
# nu=0.2;
# E = G*(2*(1+nu)) ;
# λ= E*nu/((1+nu)*(1-2*nu));
# Exx=0.2; Eyy=0; Ezz=0;
# Exy=0;   Exz=0; Eyz=0
# [ Sxx,Syy,Szz,Sxy,Sxz,Syz ] = HookesLaw3dStrain2Stress( Exx,Eyy,Ezz,Exy,Exz,Eyz,λ,G )
#
# Additional elastic constant conversions if needed:
# nu =λ/(2*(G+lamda);         Equation 8.28 Pollard       
# G = E/(2*(1+nu));           Equation 8.26 Pollard
# E = G*(2*(1+nu)) ;          Equation 8.26 Pollard (rearranged)
# λ= E*nu/((1+nu)*(1-2*nu));  Equation 8.27 Pollard
#
#  Author: Tim Davis
#  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

#Preassign
Sxx=copy(Exx)
Syy=copy(Sxx);
Szz=Ezz;
Sxy=Exy;
Sxz=Exz;
Syz=Eyz;
G2=2*G
for i=1:size(Sxx,1)
	for j=1:size(Sxx,2)
		#Equation 7.131 and 7.132 in Pollard and Fletcher 2005 Book. 
		cons=λ*(Exx[i,j]+Eyy[i,j]+Ezz[i,j]);
		Sxx[i,j]= G2*Exx[i,j]+cons; 
		Syy[i,j]= G2*Eyy[i,j]+cons;
		Szz[i,j]= G2*Ezz[i,j]+cons;
		Sxy[i,j]= G2*Exy[i,j];  
		Sxz[i,j]= G2*Exz[i,j];
		Syz[i,j]= G2*Eyz[i,j];     
	end
end
return(Sxx,Syy,Szz,Sxy,Sxz,Syz) 

end



#Multiple dispatch - working within structures
function HookesLaw3DStrain2Stress!( StrainTensor::Strains,λ,G )

	(σxx,σyy,σzz,σxy,σxz,σyz)=HookesLaw3DStrain2Stress!( StrainTensor.εxx,StrainTensor.εyy,StrainTensor.εzz,StrainTensor.εxy,StrainTensor.εxz,StrainTensor.εyz,λ,G )
	StressTensor=Stresses(σxx,σyy,σzz,σxy,σxz,σyz);
	return(StressTensor) 

end
