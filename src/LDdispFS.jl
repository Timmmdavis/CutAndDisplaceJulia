function LDdispFS(x,y,xe,ye,a,Beta,Ds,Dn,nu)
# LDdispFS: LineDisplacementInducedDisplacementFullSpace Computes the influence
#               (stress) of a single planar line crack shearing and/or
#               opening on the surronding points
#
# usage #1:
#[Stress] = LDdispFS(x,y,xe,ye,a,Beta,Ds,Dn,nu)
#
# Arguments: (input)
#      x & y  - The observation points locations in real Cartesian coords
#
#     xe & ye - The element midpoint location in real coords
#
#       a     - The elements half length
#
#       Beta  - The angle of the element away from the x-axis (radians)
#               When the normal in the -y axis down this is 0 In degrees
#               when the normal points east this is 90, west -90 and north
#               180
#
#     Dn,Ds   - The defined displacement of each element(normal and shear)
#               Dn+ is opening, Ds+ is left lateral shearing
#
#       nu    - The Poisson's ratio
#
# Arguments: (output)
# Sxx,Syy,Sxy - Is the stress caused by the movement of the dislocation at the
#               observatation points [Sxx,Syy,Sxy]
#
# Example usage:
#
#	x = [-2:0.05:2;];
#	y = [-4:0.05:0;];
#	x,y = MyModule.meshgrid(x,y);
#	dimx,dimy = size(x);
#	x=reshape(x,length(x),1);
#	y=reshape(y,length(y),1);
#   Ds=0; Dn=1;
#	(Ux,Uy)=MyModule.LDdispFS(x,y,0,-2,1,0,Ds,Dn,0.25);
#   using PyPlot
#	quiver(x, y, Ux, Uy)
#
#  Author: Tim Davis/Steve Martel
#  Copyright 2018, Tim Davis, Potsdam University
#  Modified form of Steve Martels BEM scripts from:
#  http://wwwsoesthawaiiedu/martel/MartelBEM_dir/
#  Also includes some notation from Ritz & Pollard 2012


# Define material constant used in calculating influence coefficients
con=1/(4*pi*(1-nu));
Dxb = Ds; Dyb = -Dn;
sb = sin(Beta); cb = cos(Beta);
#Define material constants used in calculating displacements.
pr1 = 1-2*nu; pr2 = 2-2*nu;

#Init array (disp at each point)
Ux= Array{Float64}(undef, length(x),1);
Uy= Array{Float64}(undef, length(x),1);

for j = 1:length(x); 
	

    # Define array of local coordinates for the observation grid relative to
    #   the midpoint and orientation of the ith element
    # Refer to (Figure 56, C&S, p 91) and eqs 451 of C&S, p 57
    XB = (x[j]-xe)*cb + (y[j]-ye)*sb;
    YB = -(x[j]-xe)*sb + (y[j]-ye)*cb;

    # Calculate derivatives of the function f(x,y), eq 525 of C&S, p 81
    #   which are used to calculate the displacement and stress components
    # It is understood that X and Y refer to XB and YB
    # First abbreviate repeated terms in the derivatives of f(x,y):
    Y2 = YB^2;
    XMa = XB-a; XPa = XB+a;
    XMa2 = XMa^2; XPa2 = XPa^2;
    R1S = XMa2 + Y2; 
    R2S = XPa2 + Y2;
	#The following derivatives are eqs. 4.5.5a thru d of C&S, p. 58.
	FF2 = con*(log(sqrt(R1S)) - log(sqrt(R2S)));	
	
	
	#Steve Martels Solution to elements lying on same plane
	#FB3 = 0 for pts colinear with element, CON*pi for pts. on element 
	#FB3 = difference of arc tangents for all other pts.
	if YB==0; 
		if abs(XB) < a;
			FF3=pi;
		else
			FF3=0;
		end
	else
		FF3 = atan(YB,XMa) - atan(YB,XPa);
	end
	FF3 = -con.*(FF3); 
    FF4 = con*(YB/R1S - YB/R2S);
    FF5 = con*(XMa/R1S - XPa/R2S);

	#Calculate the displacement components using eqs. 5.5.4 of C&S, p. 91.
	Ux[j,1] = Dxb*(-pr1*sb*FF2 + pr2*cb*FF3 + YB.*(sb*FF4 - cb*FF5))+
			  Dyb*(-pr1*cb*FF2 - pr2*sb*FF3 - YB.*(cb*FF4 + sb*FF5));
	Uy[j,1] = Dxb*(+pr1*cb*FF2 + pr2*sb*FF3 - YB.*(cb*FF4 + sb*FF5))+
			  Dyb*(-pr1*sb*FF2 + pr2*cb*FF3 - YB.*(sb*FF4 - cb*FF5));
    
end
	
return(Ux,Uy)

end
