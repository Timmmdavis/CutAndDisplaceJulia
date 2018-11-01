function LDHSInfMat(x,y,xe,ye,a,Beta,nu,Mu,DispFlag)
# LDHSInfMat: LineDisplacementStressHalfSpace Computes the influence (stress) of a single 
#			  planar line crack shearing and/or
#             opening on the surronding points (half space)
#
# usage #1:
#Sxx,Syy,Sxy,Ds,Dn = LDHSInfMat(x,y,xe,ye,a,Beta,nu,Mu,1)
# usage #2:
#Sxx,Syy,Sxy = LDHSInfMat(x,y,xe,ye,a,Beta,nu,Mu,0)
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
#x = [-2:0.5:2;];
#y = [-4:0.5:0;];
#x,y=MyModule.meshgrid(x,y);
#dimx,dimy = size(x);
#x=reshape(x,length(x),1);
#y=reshape(y,length(y),1);
#DispFlag=0;
#(SxxDs,SyyDs,SxyDs,SxxDn,SyyDn,SxyDn)=MyModule.LDHSInfMat(x,y,x,y,ones(size(x))*0.01,zeros(size(x)),0.25,1,DispFlag);
#
#  Author: Tim Davis/Steve Martel
#  Copyright 2018, Tim Davis, Potsdam University
#  Modified form of Steve Martels BEM scripts from:
#  http://wwwsoesthawaiiedu/martel/MartelBEM_dir/
#  Also includes some notation from Ritz & Pollard 2012


# Define material constant used in calculating influence coefficients
con=1/(4*pi*(1-nu));
cons=2*Mu;
Ds = 1; Dn = -1;


if DispFlag==1
	#Define material constants used in calculating displacements.
	pr1 = 1-2*nu; pr2 = 2-2*nu; pr3=3-4*nu;
	#Init array (stress at each point)
	UxDs= Array{Float64}(undef, length(x),1);
	UyDs= Array{Float64}(undef, length(x),1);
	UxDn= Array{Float64}(undef, length(x),1);
	UyDn= Array{Float64}(undef, length(x),1);
end

#Init array (stress at each point)
SxxDn= Array{Float64}(undef, length(x),length(xe));
SyyDn= Array{Float64}(undef, length(x),length(xe));
SxyDn= Array{Float64}(undef, length(x),length(xe));

SxxDs= Array{Float64}(undef, length(x),length(xe));
SyyDs= Array{Float64}(undef, length(x),length(xe));
SxyDs= Array{Float64}(undef, length(x),length(xe));

#Finding els 'j' influence on MidPoints of el 'k'
for j = 1:length(x); 
	
	sb = sin(Beta[j]); cb = cos(Beta[j]);
	s2b = sin(2*Beta[j]); c2b = cos(2*Beta[j]);
	s3b = sin(3*Beta[j]); c3b = cos(3*Beta[j]);
	s4b = sin(4*Beta[j]); c4b = cos(4*Beta[j]);

	for k=1:length(xe); #Threads.@threads #For faster runtimes...

		if (y[j]>0)
			error("Half-space solution: Z coordinates must be negative!")
		end
		
		
		# Define array of local coordinates for the observation grid relative to
		#   the midpoint and orientation of the ith element
		# Refer to (Figure 56, C&S, p 91) and eqs 451 of C&S, p 57
		XB = (x[j]-xe[k])*cb + (y[j]-ye[k])*sb;
		YB = -(x[j]-xe[k])*sb + (y[j]-ye[k])*cb;

		# Coordinates of the image dislocation
		XBi = (x[j]-xe[k])*cb - (y[j]+ye[k])*sb;		#equation 746 C&S
		YBi = (x[j]-xe[k])*sb + (y[j]+ye[k])*cb;

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
		XMa = XB-a[j]; XPa = XB+a[j];
		XMa2 = XMa^2; XPa2 = XPa^2;
		R1S = XMa2 + Y2; R1S2 = R1S^2;
		R2S = XPa2 + Y2; R2S2 = R2S^2;

		# Same thing for the image dislocation
		Y2i = YBi^2;
		XMai = XBi-a[j]; XPai = XBi+a[j];
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
		FF8i =(YBi*(1/((a[j] + XBi)^2 + YBi^2)^2 - 1/(YBi^2 + (a[j] - XBi)^2)^2 + (2*(a[j] - XBi)*(2*a[j] - 2*XBi))/(YBi^2 + (a[j] - XBi)^2)^3 - (2*(a[j] + XBi)*(2*a[j] + 2*XBi))/((a[j] + XBi)^2 + YBi^2)^3))/(2*pi*(nu - 1));

		FF9i =((a[j] - XBi)/(YBi^2 + (a[j] - XBi)^2)^2 + (a[j] + XBi)/((a[j] + XBi)^2 + YBi^2)^2)/(2*pi*(nu - 1)) - (YBi*((4*YBi*(a[j] + XBi))/((a[j] + XBi)^2 + YBi^2)^3 + (4*YBi*(a[j] - XBi))/(YBi^2 + (a[j] - XBi)^2)^3))/(2*pi*(nu - 1));


		# Calculate the stress components using eqs 555 of C&S, p 92
		SxxDsJ = cons*Ds*(2*(cb*cb)*FF4 + s2b*FF5 + YB*(c2b*FF6-s2b*FF7));
		SxxDnJ = cons*Dn*(-FF5 + YB*(s2b*FF6 + c2b*FF7));
		SyyDsJ = cons*Ds*(2*(sb*sb)*FF4 - s2b*FF5 - YB*(c2b*FF6-s2b*FF7));
		SyyDnJ = cons*Dn*(-FF5 - YB*(s2b*FF6 + c2b*FF7));
		SxyDsJ = cons*Ds*(s2b*FF4 - c2b*FF5 + YB*(s2b*FF6+c2b*FF7));
		SxyDnJ = cons*Dn*(-YB*(c2b*FF6 - s2b*FF7));

		#  Calculate IMAGE AND SUPPLEMENTAL STRESS components due to unit SHEAR and
		#  NORMAL displacement discontinuity
		SxxiDs = cons*Ds*(FF4i - 3*(c2b*FF4i - s2b*FF5i) +
		(2*y[j]*(cb - 3*c3b) + 3*YB*c2b)*FF6i +
		(2*y[j]*(sb - 3*s3b) + 3*YB*s2b)*FF7i -
		2*y[j]*(y[j]*c4b - YB*c3b)*FF8i -
		2*y[j]*(y[j]*s4b - YB*s3b)*FF9i);

		SxxiDn = cons*Dn*(FF5i + (2*y[j]*(sb - 2*s3b) +
		3*YB*s2b)*FF6i - (2*y[j]*(cb - 2*c3b) +
		3*YB*c2b)*FF7i - 2*y[j]*(y[j]*s4b - YB*s3b)*FF8i +
		2*y[j]*(y[j]*c4b - YB*c3b)*FF9i);

		SyyiDs = cons*Ds*(FF4i - (c2b*FF4i - s2b*FF5i) -
		(4*y[j]*sb*s2b - YB*c2b)*FF6i +
		(4*y[j]*sb*c2b + YB*s2b)*FF7i +
		2*y[j]*(y[j]*c4b - YB*c3b)*FF8i +
		2*y[j]*(y[j]*s4b - YB*s3b)*FF9i);

		SyyiDn = cons*Dn*(FF5i - (2*y[j]*sb - YB*s2b)*FF6i +
		(2*y[j]*cb - YB*c2b)*FF7i +
		2*y[j]*(y[j]*s4b - YB*s3b)*FF8i -
		2*y[j]*(y[j]*c4b - YB*c3b)*FF9i);

		SxyiDs = cons*Ds*(s2b*FF4i + c2b*FF5i +
		(2*y[j]*sb*(1+4*c2b) - YB*s2b)*FF6i +
		(2*y[j]*cb*(3-4*c2b) + YB*c2b)*FF7i +
		2*y[j]*(y[j]*s4b - YB*s3b)*FF8i -
		2*y[j]*(y[j]*c4b - YB*c3b)*FF9i);

		SxyiDn = cons*Dn*((4*y[j]*sb*s2b + YB*c2b)*FF6i -
		(4*y[j]*sb*c2b - YB*s2b)*FF7i -
		2*y[j]*(y[j]*c4b - YB*c3b)*FF8i -
		2*y[j]*(y[j]*s4b - YB*s3b)*FF9i);


		SxxDs[j,k]=SxxDsJ+SxxiDs;
		SxxDn[j,k]=SxxDnJ+SxxiDn;
		SyyDs[j,k]=SyyDsJ+SyyiDs;
		SyyDn[j,k]=SyyDnJ+SyyiDn;
		SxyDs[j,k]=SxyDsJ+SxyiDs;
		SxyDn[j,k]=SxyDnJ+SxyiDn;
		
		if DispFlag==1
			#The following derivatives are eqs. 4.5.5a thru d of C&S, p. 58.
			FF2 = con*(log(sqrt(R1S)) - log(sqrt(R2S)));	
			# Calculate intermediate functions Fnifor the image dislocation
			FF2i = con*(log(sqrt(R1Si)) - log(sqrt(R2Si)) ); #Equations 4.5.5 C&S	
			
			#Steve Martels Solution to elements lying on same plane
			#FB3 = 0 for pts colinear with element, CON*pi for pts. on element 
			#FB3 = difference of arc tangents for all other pts.
			if YB==0; 
				if abs(XB) < a[j];
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
			

			UxDs[j,1]=UxJ+Uxi_s;
			UxDn[j,1]=UxJ+Uxi_n;
			UyDs[j,1]=UyJ+Uyi_s;	
			UyDn[j,1]=UyJ+Uyi_n;
			
		end		

	end
    
end
	
if DispFlag==0
	return(SxxDs,SyyDs,SxyDs,SxxDn,SyyDn,SxyDn)
else
	return(SxxDs,SyyDs,SxyDs,SxxDn,SyyDn,SxyDn,UxDs,UxDn,UyDs,UyDn)
end

end
