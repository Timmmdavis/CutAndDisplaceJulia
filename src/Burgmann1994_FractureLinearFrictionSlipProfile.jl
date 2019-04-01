function Burgmann1994_FractureLinearFrictionSlipProfile(G,ν,Ts,x,a)
# Burgmann1994_FractureLinearFrictionSlipProfilet: Returns the displacement
#               of fracture walls for a fracture alligned with the x-axis 
#               that has a a linear friction profile loaded by a constant
#               shear stress. This profile is a cohesive strength that is 0
#               at the crack centre, increasing linearly to the driving
#               shear stress at the fracture tips.
#
#               See: Bürgmann, R., Pollard, D.D. and Martel, S.J., 1994.
#               Slip distributions on faults: effects of stress gradients,
#               inelastic deformation, heterogeneous host-rock stiffness,
#               and fault interaction. Journal of Structural Geology,
#               16(12), pp.1675-1690.
#
# usage #1:
# [Disp]=Burgmann1994_FractureLinearFrictionSlipProfilet(G,ν,Sxy,x)
#
# Arguments: (input)
#       ν    - Poisson's Ratio.
#
#       G    - Shear Modulus.
#
#       Ts   - The shear traction acting on the cracks faces
#
#       x     - The location of the points on the fractures walls (we
#              assume the fracture is alligned with the x-axis and extends from 1
#              to -1). 
#
#       a     - Unit half-length displacement discontinuity. This is going to have its 
#               centre at 0,0. 
#
# Arguments: (output)
#       Disp   - The Burgers vector (separation of dislocation walls) at
#               each discrete x location.
#
# Example usage:
#
#  x  = linspace(-1,1,50); #50 evenly spaced grid points 
#  G = 500; 
#  ν = 0.25;
#  Sxy= 1; 
#  a = 1; 
#  
#  [Disp]=Burgmann1994_FractureLinearFrictionSlipProfile(G,ν,Sxy,x,a);
# 
#  plot(x,Disp); xlabel('X-coord'); ylabel('Displacement of walls')
#
#  Author: Tim Davis
#  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen

x=abs.(x)
# Calculate displacement 
Disp=zeros(size(x))
for i=1:length(x)
	if a/x[i]≥1 #acosh only works above 1 in Julia, in MATLAB it returns some imaginary values
		Disp[i]=((1-ν)/G)*(2*(((Ts*(sqrt(a^2-x[i]^2))))-((Ts*(1/pi))*(((sqrt(a^2-x[i]^2))+(x[i]^2/a)*acosh(a/x[i])))))); #
	else
		Disp[i]=((1-ν)/G)*(2*(((Ts*(sqrt(a^2-x[i]^2))))-((Ts*(1/pi))*(((sqrt(a^2 -x[i]^2))+(x[i]^2/a)))))); 
	end
end
Disp=real.(Disp);

return Disp

end
