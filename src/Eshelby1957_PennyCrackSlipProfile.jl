function Eshelby1957_PennyCrackSlipProfile(G,ν,Tn,Ts,r,a)
# Eshelby1957_PennyCrackSlipProfile: Returns the displacement
#               of fracture walls for a penny shaped crack
#               loaded by a constant shear/normal traction. Crack half
#               length here is 1.
#               
#               See: Segall, P., 2010. Earthquake and volcano deformation. 
#               Princeton University Press.
#
# usage #1:
# [Disp]=Eshelby1957_PennyCrackSlipProfile(G,ν,Tn,Ts,r)
# usage #2: (Assuming this is lies in the XY plane
# [Disp]=Eshelby1957_PennyCrackSlipProfile(G,ν,Syy,Sxz/Syz,r)
#
# Arguments: (input)
#       ν    - Poisson's Ratio.
#
#       G    - Shear Modulus.
#
#       Tn    - Normal traction on the fractures face. 
#
#       Ts    - Shear traction on the fractures face. 
#
#       r     - The location of the points on the fractures walls.(radial
#              distance from centre).Assuming a=1 these should be from -1
#              to 1.
#
# 		a 	  - The crack radius
#
# Arguments: (output)
#       Disp   - The Burgers vector (separation of dislocation walls) at
#               each discrete r location.
#
# Example usage:
#
#  r  = linspace(-1,1,50); 
#  G = 500; 
#  ν = 0.25;
#  Ts= 1; 
#  Tn= 1;
#  a = 1;
# 
#  [Disp]=Eshelby1957_PennyCrackSlipProfile(G,ν,Tn,Ts,r,a)
# 
#  plot(r,Disp); xlabel('radial-coord'); ylabel('Displacement of walls') 
#
#
#  Author: Tim Davis
#  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


# Calculate displacement of penny shaped crack walls.
# These are seperated into seperate components based on related driving
# stress so the equation as a whole is easier to read.

#Shearing of crack walls # #One side of crack (eq 4.74 segall)
Ds=((4*(1-ν)*a*Ts)/(pi*(2-ν)*G)).*sqrt.(1.0 .-((r.^2)/(a^2)));
#Opening of crack walls # #Normal displacement, (eq 4.73 segall)
Dn=((2*a*(1-ν)*Tn)/(pi*G)).*sqrt.(1.0 .-((r.^2)/(a^2)));
#Both (pythag therom)
Disp=sqrt.((Ds.^2).+(Dn.^2));
# distance between relative faces;
Disp=Disp.*2;

return Disp

# ###Pressure to vol relationship for crack loaded by shear stress
# Vol=(8*Ts*a^3*(ν - 1))/(3*G*(ν - 2))
# ###Pressure to vol relationship for crack loaded by normal stress
# Vol=(2*Tn*a^3*(1 - r^2/a^2)^(1/2)*(ν - 1))/G
end