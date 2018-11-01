function LDstressFS(x,y,xe,ye,a,Beta,Ds,Dn,nu,Mu)
# LDstressFS: LineDisplacementStressFullSpace Computes the influence
#               (stress) of a single planar line crack shearing and/or
#               opening on the surronding points .
#
# usage #1:
#[Stress] = LDstressFS(x,y,xe,ye,a,Beta,Ds,Dn,nu,Mu)
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
#	(Sxx,Syy,Sxy)=MyModule.LDstressFS(x,y,0,-2,1,0,Ds,Dn,0.25,1);
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

#Init array (stress at each point)
Sxx= Array{Float64}(undef, length(x),1);
Syy= Array{Float64}(undef, length(x),1);
Sxy= Array{Float64}(undef, length(x),1);

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
    R1S = XMa2 + Y2; R1S2 = R1S^2;
    R2S = XPa2 + Y2; R2S2 = R2S^2;

    FF4 = con*(YB/R1S - YB/R2S);
    FF5 = con*(XMa/R1S - XPa/R2S);
    # The following derivatives are eqs 553a and b of C&S, p 91
    FF6 = con*((XMa2 - Y2)/R1S2 - (XPa2 - Y2)/R2S2);
    FF7 = 2*con*YB*(XMa/R1S2 - XPa/R2S2);


    # Calculate the stress components using eqs 555 of C&S, p 92
    Sxx[j,1] =  cons*Dxb*(2*(cb*cb)*FF4 + s2b*FF5 + YB*(c2b*FF6-s2b*FF7))+
    cons*Dyb*(-FF5 + YB*(s2b*FF6 + c2b*FF7));
    Syy[j,1] =  cons*Dxb*(2*(sb*sb)*FF4 - s2b*FF5 - YB*(c2b*FF6-s2b*FF7))+
    cons*Dyb*(-FF5 - YB*(s2b*FF6 + c2b*FF7));
    Sxy[j,1] =  cons*Dxb*(s2b*FF4 - c2b*FF5 + YB*(s2b*FF6+c2b*FF7))+
    cons*Dyb*(-YB*(c2b*FF6 - s2b*FF7));
    
end
	
return(Sxx,Syy,Sxy)

end
