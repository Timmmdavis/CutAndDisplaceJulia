function HookesLaw3DStrain2Stress( Exx,Eyy,Ezz,Exy,Exz,Eyz,λ,G )
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

#Equation 7.131 and 7.132 in Pollard and Fletcher 2005 Book. 
Sxx = 2*G.*Exx.+(λ.*(Exx.+Eyy.+Ezz)); 
Syy = 2*G.*Eyy.+(λ.*(Exx.+Eyy.+Ezz));
Szz = 2*G.*Ezz.+(λ.*(Exx.+Eyy.+Ezz));
Sxy = 2*G.*Exy;  
Sxz = 2*G.*Exz;
Syz = 2*G.*Eyz;     

return(Sxx,Syy,Szz,Sxy,Sxz,Syz) 

end



#Multiple dispatch - working within structures
function HookesLaw3dStrain2Stress( StrainTensor::Strains,λ,G )

	(σxx,σyy,σzz,σxy,σxz,σyz)=HookesLaw3dStrain2Stress( StrainTensor.εxx,StrainTensor.εyy,StrainTensor.εzz,StrainTensor.εxy,StrainTensor.εxz,StrainTensor.εyz,λ,G )
	StressTensor=Stresses(σxx,σyy,σzz,σxy,σxz,σyz);
	return(StressTensor) 

end
