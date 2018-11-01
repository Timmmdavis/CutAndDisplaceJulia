function LDdispHS(x,y,xe,ye,a,Beta,Ds,Dn,nu)
# LDdispFS: LineDisplacementInducedDisplacementFullSpace Computes the influence
#               (stress) of a single planar line crack shearing and/or
#               opening on the surronding points
#
# usage #1:
#[Stress] = LDdispHS(x,y,xe,ye,a,Beta,Ds,Dn,nu)
#
# Arguments: (input)
#      x & y  - The observation points locations in real Cartesian coords
#              (y must be below 0)
#
#     xe & ye - The element midpoint location in real coords
#              (must have a ye below 0)
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
#	(Ux,Uy)=MyModule.LDdispHS(x,y,0,-2,1,0,Ds,Dn,0.25);
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
s2b = sin(2*Beta); c2b = cos(2*Beta);
s3b = sin(3*Beta); c3b = cos(3*Beta);

#Define material constants used in calculating displacements.
pr1 = 1-2*nu; pr2 = 2-2*nu; pr3=3-4*nu;

#Init array (disp at each point)
Ux= Array{Float64}(undef, length(x),1);
Uy= Array{Float64}(undef, length(x),1);

for j = 1:length(x); 
	

    # Define array of local coordinates for the observation grid relative to
    #   the midpoint and orientation of the ith element
    # Refer to (Figure 56, C&S, p 91) and eqs 451 of C&S, p 57
    XB = (x[j]-xe)*cb + (y[j]-ye)*sb;
    YB = -(x[j]-xe)*sb + (y[j]-ye)*cb;
	
    # Coordinates of the image dislocation
    XBi = (x[j]-xe)*cb - (y[j]+ye)*sb;		#equation 746 C&S
    YBi = (x[j]-xe)*sb + (y[j]+ye)*cb;
	
    # Fix roundoff errors in Ybi and Yb from trig function problems
    for i=1:length(YBi)
        if abs(YBi)<1e-10
            YBi=0;
        end
    end
    for i=1:length(YB)
        if abs(YB)<1e-10
            YB=0;
        end
    end	

    # Calculate derivatives of the function f(x,y), eq 525 of C&S, p 81
    #   which are used to calculate the displacement and stress components
    # It is understood that X and Y refer to XB and YB
    # First abbreviate repeated terms in the derivatives of f(x,y):
    Y2 = YB^2;
    XMa = XB-a; XPa = XB+a;
    XMa2 = XMa^2; XPa2 = XPa^2;
    R1S = XMa2 + Y2; 
    R2S = XPa2 + Y2;
	
    # Same thing for the image dislocation
    Y2i = YBi^2;
    XMai = XBi-a; XPai = XBi+a;
    XMa2i = XMai^2; XPa2i = XPai^2;
    R1Si = XMa2i + Y2i; R1S2i = R1Si^2;
    R2Si = XPa2i + Y2i; R2S2i = R2Si^2;	
	
	#The following derivatives are eqs. 4.5.5a thru d of C&S, p. 58.
	FF2 = con*(log(sqrt(R1S)) - log(sqrt(R2S)));	
	# Calculate intermediate functions Fnifor the image dislocation
	FF2i = con*(log(sqrt(R1Si)) - log(sqrt(R2Si)) ); #Equations 4.5.5 C&S
	
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
	#And same for image dislocation
	if YBi==0; 
		if abs(XBi) < a;
			FF3i=pi;
		else
			FF3i=0;
		end
	else
		FF3i = atan(YBi,XMai) - atan(YBi,XPai);
	end
	FF3i = -con.*(FF3i); 
	
	
    FF4 = con*(YB/R1S - YB/R2S);
    FF5 = con*(XMa/R1S - XPa/R2S);
	
    FF4i = con*(YBi/R1Si - YBi/R2Si);
    FF5i = con*(XMai/R1Si - XPai/R2Si);	
	
	
	# The halfspace examples of eqs 553a and b of C&S, p 91
    # See Appendix A of: Martel, SJ and Langley, JS, 2006 Propagation of
    # normal faults to the surface in basalt, Koae fault system, Hawaii
    # Journal of Structural Geology, 28(12), pp2123-2143
    FF6i = con*((XMa2i - Y2i)/R1S2i - (XPa2i - Y2i)/R2S2i);
    FF7i = 2*con*YBi*(XMai/R1S2i - XPai/R2S2i);

	#Calculate the displacement components using eqs. 5.5.4 of C&S, p. 91.
	UxJ = Dxb*(-pr1*sb*FF2 + pr2*cb*FF3 + YB.*(sb*FF4 - cb*FF5))+
		  Dyb*(-pr1*cb*FF2 - pr2*sb*FF3 - YB.*(cb*FF4 + sb*FF5));
	UyJ = Dxb*(+pr1*cb*FF2 + pr2*sb*FF3 - YB.*(cb*FF4 + sb*FF5))+
		  Dyb*(-pr1*sb*FF2 + pr2*cb*FF3 - YB.*(sb*FF4 - cb*FF5));
    
	#  See equations 7.4.8 and 7.4.9 in Crouch and Starfield
	#  Calculate image and supplemental displacement components due to unit shear displacement discontinuity
	Uxi_s =Dxb* (pr1*sb*FF2i - pr2*cb*FF3i +
    (pr3*(y[j]*s2b - YB*sb) + 2*y[j]*s2b)*FF4i +
    (pr3*(y[j]*c2b - YB*cb) - y[j]*(1-2*c2b))*FF5i +
    2*y[j]*(y[j]*s3b - YB*s2b)*FF6i -
    2*y[j]*(y[j]*c3b - YB*c2b)*FF7i);

	Uyi_s =Dxb* (-pr1*cb*FF2i - pr2*sb*FF3i-
    (pr3*(y[j]*c2b - YB*cb) + y[j]*(1-2*c2b))*FF4i+  
    (pr3*(y[j]*s2b - YB*sb) - 2*y[j]*s2b)*FF5i+	 
    2*y[j]*(y[j]*c3b - YB*c2b)*FF6i+
    2*y[j]*(y[j]*s3b - YB*s2b)*FF7i);     

	#Calculate image and supplemental displacement components due to unit normal displacement discontinuity
	Uxi_n =Dyb*(pr1*cb*FF2i + pr2*sb*FF3i-
    (pr3*(y[j]*c2b - YB*cb) - y[j])*FF4i+
    pr3*(y[j]*s2b - YB*sb)*FF5i-
    2*y[j]*(y[j]*c3b - YB*c2b)*FF6i-
    2*y[j]*(y[j]*s3b - YB*s2b)*FF7i);    

	Uyi_n =Dyb* (pr1*sb*FF2i - pr2*cb*FF3i-
    pr3*(y[j]*s2b - YB*sb)*FF4i-
    (pr3*(y[j]*c2b - YB*cb) + y[j])*FF5i+
    2*y[j]*(y[j]*s3b - YB*s2b)*FF6i-
    2*y[j]*(y[j]*c3b - YB*c2b)*FF7i);

	Uxi=Uxi_s+Uxi_n;
	Uyi=Uyi_s+Uyi_n;

	Ux[j,1]=UxJ+Uxi;
	Uy[j,1]=UyJ+Uyi;	
	
end
	
return(Ux,Uy)

end
