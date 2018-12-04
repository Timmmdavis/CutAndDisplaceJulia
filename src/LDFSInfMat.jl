function LDFSInfMat(x,y,xe,ye,a,Beta,nu,Mu,DispFlag)
# LDFSInfMat: LineDisplacementStressFullSpace Creates inf matricies
# 			  faster as most constants are the same for Ds and Dn. 
#			  Inf matricies are one unit of each. 
#
# usage #1:
#	Sxx,Syy,Sxy,Ds,Dn = LDFSInfMat(x,y,xe,ye,a,Beta,nu,Mu,1)
# usage #2:
#	Sxx,Syy,Sxy = LDFSInfMat(x,y,xe,ye,a,Beta,nu,Mu,0)
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
#       nu    - The Poisson's ratio
#
#       Mu     - The shear modulus
#
#     DispFlag - Create displacement inf mats also
#
# Arguments: (output)
# Sxx,Syy,Sxy - Is the stress caused by the movement of the dislocation at the
#               observatation points [Sxx,Syy,Sxy]. Dn means those caused by a
#				magnitude 1 of opening. Ds for shearing
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
	pr1 = 1-2*nu; pr2 = 2-2*nu;
	#Init array (stress at each point)
	UxDs= Array{Float64}(undef, length(x),length(x));
	UyDs= Array{Float64}(undef, length(x),length(x));
	UxDn= Array{Float64}(undef, length(x),length(x));
	UyDn= Array{Float64}(undef, length(x),length(x));
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

	for k=1:length(xe); #Threads.@threads #For faster runtimes...
		
		# Define array of local coordinates for the observation grid relative to
		#   the midpoint and orientation of the ith element
		# Refer to (Figure 56, C&S, p 91) and eqs 451 of C&S, p 57
		XB = (x[j]-xe[k])*cb + (y[j]-ye[k])*sb;
		YB = -(x[j]-xe[k])*sb + (y[j]-ye[k])*cb;

		# Calculate derivatives of the function f(x,y), eq 525 of C&S, p 81
		#   which are used to calculate the displacement and stress components
		# It is understood that X and Y refer to XB and YB
		# First abbreviate repeated terms in the derivatives of f(x,y):
		Y2 = YB^2;
		XMa = XB-a[j]; XPa = XB+a[j];
		XMa2 = XMa^2; XPa2 = XPa^2;
		R1S = XMa2 + Y2; R1S2 = R1S^2;
		R2S = XPa2 + Y2; R2S2 = R2S^2;

		FF4 = con*(YB/R1S - YB/R2S);
		FF5 = con*(XMa/R1S - XPa/R2S);
		# The following derivatives are eqs 553a and b of C&S, p 91
		FF6 = con*((XMa2 - Y2)/R1S2 - (XPa2 - Y2)/R2S2);
		FF7 = 2*con*YB*(XMa/R1S2 - XPa/R2S2);

		# Calculate the stress components using eqs 555 of C&S, p 92
		SxxDs[j,k] = cons*Ds*(2*(cb*cb)*FF4 + s2b*FF5 + YB*(c2b*FF6-s2b*FF7));
		SxxDn[j,k] = cons*Dn*(-FF5 + YB*(s2b*FF6 + c2b*FF7));
		SyyDs[j,k] = cons*Ds*(2*(sb*sb)*FF4 - s2b*FF5 - YB*(c2b*FF6-s2b*FF7));
		SyyDn[j,k] = cons*Dn*(-FF5 - YB*(s2b*FF6 + c2b*FF7));
		SxyDs[j,k] = cons*Ds*(s2b*FF4 - c2b*FF5 + YB*(s2b*FF6+c2b*FF7));
		SxyDn[j,k] = cons*Dn*(-YB*(c2b*FF6 - s2b*FF7));
		
		if DispFlag==1
			#The following derivatives are eqs. 4.5.5a thru d of C&S, p. 58.
			FF2 = con*(log(sqrt(R1S)) - log(sqrt(R2S)));	
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
			#Calculate the displacement components using eqs. 5.5.4 of C&S, p. 91.
			UxDs[j,k] = Ds*(-pr1*sb*FF2 + pr2*cb*FF3 + YB.*(sb*FF4 - cb*FF5));
			UxDn[j,k] = Dn*(-pr1*cb*FF2 - pr2*sb*FF3 - YB.*(cb*FF4 + sb*FF5));
			UyDs[j,k] = Ds*(+pr1*cb*FF2 + pr2*sb*FF3 - YB.*(cb*FF4 + sb*FF5));
			UyDs[j,k] = Dn*(-pr1*sb*FF2 + pr2*cb*FF3 - YB.*(sb*FF4 - cb*FF5));	
		end
    
	end
	
end
	
if DispFlag==0
	return(SxxDs,SyyDs,SxyDs,SxxDn,SyyDn,SxyDn)
else
	return(SxxDs,SyyDs,SxyDs,SxxDn,SyyDn,SxyDn,UxDs,UxDn,UyDs,UyDn)
end

end
