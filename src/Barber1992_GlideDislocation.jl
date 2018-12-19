function Barber1992_GlideDislocation(k,mu,X,Y,a,b,nu)
# Barber1992_GlideDislocation: Returns Cartesian displacements and stresses
#              at grid points for a displacement discontinuity with a unit
#              Burger's vector (B=1) and unit half length (a=1). The
#              displacement discontinuity extends along the x-axis from 
#              x = -a to x = +a. 
#              Note a positive Burgers vector here drives a right lateral
#              shearing. 
#
#              The solution is based on equations for a glide
#              dislocation from Barber: Elasticity (1992) and Pollard and
#              Fletcher (2005). 
#
#
#
# usage #1: [X,Y,Sxx,Syy,Sxy,Ux,Uy] =
# Barber1992_GlideDislocation(k,mu,X,Y,a,b,nu);
#
# Arguments: (input)
#       k     - Kolosov's constant for plane strain.
#
#       mu    - Shear Modulus.
#
#  minx,maxx  - The bounds of the grid that stresses and displacements are
#              found on. (Y and X).
#
#  spacing    - The spacing between the points on this grid.
#
#       a     - The half length of the dislocation.
#
#       b     - The Burgers vector (separation of dislocation walls).
#
#       nu    - Poisson's Ratio.
#
#
# Arguments: (output)
#
#       X,Y    - X and Y locations of the observation points.
#
#  Sxx,Syy,Sxy - 2D stress tensor components returned on a grid.
#
#       Ux,Uy  - Displacement of stress tensor components returned on a
#               grid.
#
#
# Example usage:
#
#  nu = 0.25;
#  k = 3-4*nu; 
#  mu = 500;
#  spacing=0.1;
#  minx=-4; maxx=4;
#  [X,Y] = meshgrid(minx:spacing:maxx);
#  a = 1;  
#  b=0.0001; 
#  
#  [X,Y,Sxx,Syy,Sxy,Ux,Uy] =...
#   Barber1992_GlideDislocation(k,mu,X,Y,a,b,nu);
# 
#  quiver(X(:),Y(:),Ux(:),Uy(:))
#  DrawContourFPlots2d( X,Y,[], Sxx,Syy,Sxy )
#
#
#  Author: Tim Davis
#  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen
#  Modified from Steve Martel's fracture mechanics homework

Sxx= zeros(length(X),1);
Syy= zeros(length(X),1);
Sxy= zeros(length(X),1);
Ux= zeros(length(X),1);
Uy= zeros(length(X),1);

for i=1:length(X)
	#Define larger grid that covers where both dislocations will be, this is
	#used to calculate stress once. 
	(Sxxp,Syyp,Sxyp,Uxp,Uyp)=CompD(mu,b,k,nu,X[i]-a,Y[i]);
	(Sxxn,Syyn,Sxyn,Uxn,Uyn)=CompD(mu,b,k,nu,X[i]+a,Y[i]);

	Ux[i]=Uxp-Uxn;
	Uy[i]=Uyp-Uyn;

	Sxx[i]=Sxxp-Sxxn;
	Syy[i]=Syyp-Syyn;
	Sxy[i]=Sxyp-Sxyn;
end


return(X,Y,Sxx,Syy,Sxy,Ux,Uy)
end

function CompD(mu,b,k,nu,X,Y)

#Set up vars
r = sqrt(X^2 + Y^2);
sint = Y/r; 
cost = X/r;  
theta = atan(Y,X);
    
    
# Calculate polar rt stress components due to unit displacement discontinuity
# Positive end Negative end
# Barber page 201 equations 1324-1325
Srr = ((2*mu)*b*sint)/(pi*(k+1)*r);
Stt = Srr; 
Srt = -((2*mu)*b*cost)/(pi*(k+1)*r);

#Converting the components into Cartesian tensor
#Calling external function
(Sxx,Syy,Sxy)=TensorTransformation2D(Srr,Stt,Srt,cos(theta),cos((pi/2)-theta) )


#Displacement equations, see Pollard and Fletcher 2005 Eq 8.36-8.37. 
#Equations from Fig 08_10
lambda = (2*mu*nu)/(1-2*nu); 
#Constants
c1 = (0.5*mu)/(lambda+2*mu); 
c2 = (lambda+mu)/(lambda+2*mu);
#Equations
#Ux
Ux1 = -(b/(2*pi))*atan.(Y,X);
Ux2 = -(b/(2*pi))*c2*(X*Y)/(X^2 +Y^2);

Ux=Ux1+Ux2; bbb=abs(Ux)==Inf; 
if bbb==true
	Ux=0;
end

#Uy
Uy1 = -(b/(2*pi))*(-c1)*log.(X^2+Y^2);
Uy2 = -(b/(2*pi))*c2*(Y^2)/(X^2+Y^2);

Uy=Uy1+Uy2; bbb=abs(Uy)==Inf; 
if bbb==true
	Uy=0;
end
return(Sxx,Syy,Sxy,Ux,Uy)
end
