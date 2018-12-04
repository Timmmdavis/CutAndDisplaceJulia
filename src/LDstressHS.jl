function LDstressHS(x,y,xe,ye,a,Beta,Ds,Dn,nu,Mu)
# LDstressHS: LineDisplacementStressHalfSpace Computes the influence (stress) of a single 
#			  planar line crack shearing and/or
#             opening on the surronding points (half space)
#
# usage #1:
#[Stress] = LDstressHS(x,y,xe,ye,a,Beta,Ds,Dn,nu,Mu)
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
#       Mu     - The shear modulus
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
#	(Sxx,Syy,Sxy)=MyModule.LDstressHS(x,y,0,-2,1,0,Ds,Dn,0.25,1);
#   using PyPlot
#	scatter(x, y,15,Sxx);
#	cbar = colorbar()
#
#  Author: Tim Davis/Steve Martel
#  Copyright 2018, Tim Davis, Potsdam University
#  Modified form of Steve Martels BEM scripts from:
#  http://wwwsoesthawaiiedu/martel/MartelBEM_dir/
#  Also includes some notation from Ritz & Pollard 2012


# Define material constant used in calculating influence coefficients
con=1/(4*pi*(1-nu));
cons=2*Mu;
Dxb = Ds; Dyb = -Dn;
sb = sin(Beta); cb = cos(Beta);
s2b = sin(2*Beta); c2b = cos(2*Beta);
s3b = sin(3*Beta); c3b = cos(3*Beta);
s4b = sin(4*Beta); c4b = cos(4*Beta);

#Init array (stress at each point)
Sxx= Array{Float64}(undef, length(x),1);
Syy= Array{Float64}(undef, length(x),1);
Sxy= Array{Float64}(undef, length(x),1);

for j = 1:length(x); 
	
	
	
    if (y[j]>0)
        error("Half-space solution: Z coordinates must be negative!")
    end
	
	
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
    R1S = XMa2 + Y2; R1S2 = R1S^2;
    R2S = XPa2 + Y2; R2S2 = R2S^2;

    # Same thing for the image dislocation
    Y2i = YBi^2;
    XMai = XBi-a; XPai = XBi+a;
    XMa2i = XMai^2; XPa2i = XPai^2;
    R1Si = XMa2i + Y2i; R1S2i = R1Si^2;
    R2Si = XPa2i + Y2i; R2S2i = R2Si^2;

    FF4 = con*(YB/R1S - YB/R2S);
    FF5 = con*(XMa/R1S - XPa/R2S);
    # The following derivatives are eqs 553a and b of C&S, p 91
    FF6 = con*((XMa2 - Y2)/R1S2 - (XPa2 - Y2)/R2S2);
    FF7 = 2*con*YB*(XMa/R1S2 - XPa/R2S2);


    FF4i = con*(YBi/R1Si - YBi/R2Si);
    FF5i = con*(XMai/R1Si - XPai/R2Si);


    # The halfspace examples of eqs 553a and b of C&S, p 91
    # See Appendix A of: Martel, SJ and Langley, JS, 2006 Propagation of
    # normal faults to the surface in basalt, Koae fault system, Hawaii
    # Journal of Structural Geology, 28(12), pp2123-2143
    FF6i = con*((XMa2i - Y2i)/R1S2i - (XPa2i - Y2i)/R2S2i);
    FF7i = 2*con*YBi*(XMai/R1S2i - XPai/R2S2i);

    #*Tim* I used MATLABs symbolic to find these not eq's A3 and A4 of Martel
    # Used EqA1 on variable FF7i (expanded)
    FF8i =(YBi*(1/((a + XBi)^2 + YBi^2)^2 - 1/(YBi^2 + (a - XBi)^2)^2 + (2*(a - XBi)*(2*a - 2*XBi))/(YBi^2 + (a - XBi)^2)^3 - (2*(a + XBi)*(2*a + 2*XBi))/((a + XBi)^2 + YBi^2)^3))/(2*pi*(nu - 1));

    FF9i =((a - XBi)/(YBi^2 + (a - XBi)^2)^2 + (a + XBi)/((a + XBi)^2 + YBi^2)^2)/(2*pi*(nu - 1)) - (YBi*((4*YBi*(a + XBi))/((a + XBi)^2 + YBi^2)^3 + (4*YBi*(a - XBi))/(YBi^2 + (a - XBi)^2)^3))/(2*pi*(nu - 1));


    # Calculate the stress components using eqs 555 of C&S, p 92
    SxxJ =  cons*Dxb*(2*(cb*cb)*FF4 + s2b*FF5 + YB*(c2b*FF6-s2b*FF7))+
    cons*Dyb*(-FF5 + YB*(s2b*FF6 + c2b*FF7));
    SyyJ =  cons*Dxb*(2*(sb*sb)*FF4 - s2b*FF5 - YB*(c2b*FF6-s2b*FF7))+
    cons*Dyb*(-FF5 - YB*(s2b*FF6 + c2b*FF7));
    SxyJ =  cons*Dxb*(s2b*FF4 - c2b*FF5 + YB*(s2b*FF6+c2b*FF7))+
    cons*Dyb*(-YB*(c2b*FF6 - s2b*FF7));



    #  Calculate IMAGE AND SUPPLEMENTAL STRESS components due to unit SHEAR and
    #  NORMAL displacement discontinuity
    Sxxi_s = cons*Dxb*(FF4i - 3*(c2b*FF4i - s2b*FF5i) +
    (2*y[j]*(cb - 3*c3b) + 3*YB*c2b)*FF6i +
    (2*y[j]*(sb - 3*s3b) + 3*YB*s2b)*FF7i -
    2*y[j]*(y[j]*c4b - YB*c3b)*FF8i -
    2*y[j]*(y[j]*s4b - YB*s3b)*FF9i);

    Sxxi_n = cons*Dyb*(FF5i + (2*y[j]*(sb - 2*s3b) +
    3*YB*s2b)*FF6i - (2*y[j]*(cb - 2*c3b) +
    3*YB*c2b)*FF7i - 2*y[j]*(y[j]*s4b - YB*s3b)*FF8i +
    2*y[j]*(y[j]*c4b - YB*c3b)*FF9i);

    Syyi_s = cons*Dxb*(FF4i - (c2b*FF4i - s2b*FF5i) -
    (4*y[j]*sb*s2b - YB*c2b)*FF6i +
    (4*y[j]*sb*c2b + YB*s2b)*FF7i +
    2*y[j]*(y[j]*c4b - YB*c3b)*FF8i +
    2*y[j]*(y[j]*s4b - YB*s3b)*FF9i);

    Syyi_n = cons*Dyb*(FF5i - (2*y[j]*sb - YB*s2b)*FF6i +
    (2*y[j]*cb - YB*c2b)*FF7i +
    2*y[j]*(y[j]*s4b - YB*s3b)*FF8i -
    2*y[j]*(y[j]*c4b - YB*c3b)*FF9i);

    Sxyi_s = cons*Dxb*(s2b*FF4i + c2b*FF5i +
    (2*y[j]*sb*(1+4*c2b) - YB*s2b)*FF6i +
    (2*y[j]*cb*(3-4*c2b) + YB*c2b)*FF7i +
    2*y[j]*(y[j]*s4b - YB*s3b)*FF8i -
    2*y[j]*(y[j]*c4b - YB*c3b)*FF9i);

    Sxyi_n = cons*Dyb*((4*y[j]*sb*s2b + YB*c2b)*FF6i -
    (4*y[j]*sb*c2b - YB*s2b)*FF7i -
    2*y[j]*(y[j]*c4b - YB*c3b)*FF8i -
    2*y[j]*(y[j]*s4b - YB*s3b)*FF9i);


    Sxx[j,1]=SxxJ+Sxxi_s+Sxxi_n;
    Syy[j,1]=SyyJ+Syyi_s+Syyi_n;
    Sxy[j,1]=SxyJ+Sxyi_s+Sxyi_n;

    
end
	
return(Sxx,Syy,Sxy)

end
