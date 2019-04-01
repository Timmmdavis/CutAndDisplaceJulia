function PollardSegall1987_FractureSlipProfile(G,ν,σyy,σxy,x,a)
# PollardSegall1987_FractureSlipProfile: Returns the displacement
#               of fracture walls for a fracture alligned with the x-axis 
#               loaded by a constant shear/normal stress. 
#               
#               See: Pollard, D.D. and Segall, P., 1987. Theoretical
#               displacements and stresses near fractures in rock: with
#               applications to faults, joints, veins, dikes, and solution
#               surfaces. Fracture mechanics of rock, 277(349), pp.277-349.
#
# usage #1:
# [Disp]=PollardSegall1987_FractureSlipProfile(G,ν,Syy,Sxy,x)
#
# Arguments: (input)
#       ν    - Poisson's Ratio.
#
#       G    - Shear Modulus.
#
#       Syy   - The normal stress at infinity. (traction Tn)
#
#       Sxy   - The shear stress at infinity.  (traction -Ts) 
#
#       x     - The location of the points on the fractures walls (we
#              assume the fracture is alligned with the x-axis and extends from 1
#              to -1). 
#
#       a     - Unit half-length displacement discontiutity. This is going to have its 
#               centre at 0,0. 
#
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
#  Syy= 1;
#  a=1;
# 
#  [Disp]=PollardSegall1987_FractureSlipProfile(G,ν,Syy,Sxy,x,a);
# 
#  plot(x,Disp); xlabel('X-coord'); ylabel('Displacement of walls') 
#
#  Author: Tim Davis
#  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


# Calculate displacement eq 8.35a and 8.35b Pollard and Segall 1984
# These are seperated into seperate components based on related driving stress so the equation as a whole is easier to read 
Disp=(2*(1-ν)/G)*sqrt.(a.^2 .-x.^2).*(σyy+σxy);
Mode3Disp=((2*(1-ν)/G).*sqrt.(a.^2.0 .-x.^2))/1 .-ν;
return Disp,Mode3Disp

end
###Pressure vol relationship
# #Elastic constant
# Vol=-(pi*(Tn+Ts)*a^2*(ν - 1))/G