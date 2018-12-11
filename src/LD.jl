function LD(x,y,xe,ye,a,Beta,Ds,Dn,nu,Mu,DispFlag,StressFlag,HSflag)
# LDFSInfMat: Creates inf matricies (or effect of single dislocation) on stress and displacment
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
#       Mu     - The shear modulus (Only needed for displacements)
#
# 		HSflag - 		if 1 use half space formulation
#
# 		DispFlag - 		if 1 compute displacement
#
# 		StressFlag - 	if 1 compute stress
#
#
# Arguments: (output): note the size of the outputs depend on the input flags
# Sxx,Syy,Sxy - Is the stress caused by the movement of the dislocation at the
#               observatation points [Sxx,Syy,Sxy]. 
#
# Ux Uy 	  - Is the displacement caused by the movement of the dislocation at the
#               observatation points
#
#To get full mats (on output):
# Sxx=SxxDs.+SxxDn;
# Syy=SyyDs.+SyyDn;
# Sxy=SxyDs.+SxyDn;
# Ux=UxDs.+UyDs;
# Uy=UxDn.+UyDn;

# # Example usage (creating infMats):
# x = [-2:0.5:2;];
# y = [-4:0.5:0;];
# x,y=MyModule.meshgrid(x,y);
# dimx,dimy = size(x);
# x=reshape(x,length(x),1);
# y=reshape(y,length(y),1);
# xe=copy(x);
# ye=copy(y);
# DispFlag=0;
# StressFlag=1;
# HSflag=0;
# Ds=ones(size(xe));
# Dn=ones(size(xe));
#(SxxDs,SyyDs,SxyDs,SxxDn,SyyDn,SxyDn)=MyModule.LD(x,y,xe,ye,ones(size(x))*0.01,zeros(size(x)),0.25,1,Ds,Dn,DispFlag,StressFlag,HSflag);

CorrectDimsFlg= size(xe)==size(ye) 	&&
				size(xe)==size(a) 	&&
				size(xe)==size(Beta)&&
				size(xe)==size(Ds) 	&&
				size(xe)==size(Dn);
if CorrectDimsFlg!=1
	error("Element size inputs must be the same dimensions")
end
 
# Define material constant used in calculating influence coefficients
con=1/(4*pi*(1-nu));
cons=2*Mu;

Dn=0 .-Dn; #flip

if StressFlag==1
	#Init array (stress at each point)
	SxxDn= Array{Float64}(undef, length(x),length(xe));
	SyyDn= Array{Float64}(undef, length(x),length(xe));
	SxyDn= Array{Float64}(undef, length(x),length(xe));
	SxxDs= Array{Float64}(undef, length(x),length(xe));
	SyyDs= Array{Float64}(undef, length(x),length(xe));
	SxyDs= Array{Float64}(undef, length(x),length(xe));
end

if DispFlag==1
	#Define material constants used in calculating displacements.
	pr1 = 1-2*nu; pr2 = 2-2*nu; pr3=3-4*nu;
	#Init array (stress at each point)
	UxDs= Array{Float64}(undef, length(x),length(xe));
	UyDs= Array{Float64}(undef, length(x),length(xe));
	UxDn= Array{Float64}(undef, length(x),length(xe));
	UyDn= Array{Float64}(undef, length(x),length(xe));
end

for k = 1:length(xe); #for every element 

	sb = sin(Beta[k]); cb = cos(Beta[k]);
	s2b = sin(2*Beta[k]); c2b = cos(2*Beta[k]);
	
	if HSflag==1
		s3b = sin(3*Beta[k]); c3b = cos(3*Beta[k]);
		s4b = sin(4*Beta[k]); c4b = cos(4*Beta[k]);
	end
	
	for j=1:length(x); #Threads.@threads #For faster runtimes... #For every point in space

		# Define array of local coordinates for the observation grid relative to
		#   the midpoint and orientation of the ith element
		# Refer to (Figure 56, C&S, p 91) and eqs 451 of C&S, p 57
		xMxe=x[j]-xe[k];
		yMye=y[j]-ye[k];
		XB =  (xMxe)*cb + (yMye)*sb;
		YB = -(xMxe)*sb + (yMye)*cb;
		if HSflag==1
			if (y[j]>0)
				error("Half-space solution: Z coordinates must be negative!")
			end		
			yAye=y[j]+ye[k];
			# Coordinates of the image dislocation
			XBi = (xMxe)*cb - (yAye)*sb;		#equation 746 C&S
			YBi = (xMxe)*sb + (yAye)*cb;
		end
		
		if StressFlag==1
			# Calculate derivatives of the function f(x,y), eq 525 of C&S, p 81
			#   which are used to calculate the displacement and stress components
			# It is understood that X and Y refer to XB and YB
			(Y2,XMa,XPa,XMa2,XPa2,R1S,R2S,R1S2,R2S2)=AbbreviateTerms(YB,XB,a[k])
			if HSflag==1
				# Same thing for the image dislocation
				(Y2i,XMai,XPai,XMa2i,XPa2i,R1Si,R2Si,R1S2i,R2S2i)=AbbreviateTerms(YBi,XBi,a[k])
			end
		end
			
		if DispFlag==1
			(FF2,FF3)=CompFF2AndFF3(R1S,R2S,con,XB,YB,XMa,XPa,a[k])
			#And same for image dislocation
			if HSflag==1
				# Calculate intermediate functions Fnifor the image dislocation
				(FF2i,FF3i)=CompFF2AndFF3(R1Si,R2Si,con,XBi,YBi,XMai,XPai,a[k])
			end
		end		
		
		if StressFlag==1
			(FF4,FF5,FF6,FF7)=CompFF4ToFF7(R1S,R2S,YB,XMa,XPa,XMa2,XPa2,R1S2,R2S2,Y2,con)
			if HSflag==1
					(FF4i,FF5i,FF6i,FF7i)=CompFF4ToFF7(R1Si,R2Si,YBi,XMai,XPai,XMa2i,XPa2i,R1S2i,R2S2i,Y2i,con)
					
					# The halfspace examples of eqs 553a and b of C&S, p 91
					# See Appendix A of: Martel, SJ and Langley, JS, 2006 Propagation of
					# normal faults to the surface in basalt, Koae fault system, Hawaii
					# Journal of Structural Geology, 28(12), pp2123-2143

					#*Tim* I used MATLABs symbolic to find these not eq's A3 and A4 of Martel
					# Used EqA1 on variable FF7i (expanded)
					FF8i =(YBi*(1/((a[k] + XBi)^2 + YBi^2)^2 - 1/(YBi^2 + (a[k] - XBi)^2)^2 + (2*(a[k] - XBi)*(2*a[k] - 2*XBi))/(YBi^2 + (a[k] - XBi)^2)^3 - (2*(a[k] + XBi)*(2*a[k] + 2*XBi))/((a[k] + XBi)^2 + YBi^2)^3))/(2*pi*(nu - 1));
					FF9i =((a[k] - XBi)/(YBi^2 + (a[k] - XBi)^2)^2 + (a[k] + XBi)/((a[k] + XBi)^2 + YBi^2)^2)/(2*pi*(nu - 1)) - (YBi*((4*YBi*(a[k] + XBi))/((a[k] + XBi)^2 + YBi^2)^3 + (4*YBi*(a[k] - XBi))/(YBi^2 + (a[k] - XBi)^2)^3))/(2*pi*(nu - 1));
			end

			# Calculate the stress components using eqs 555 of C&S, p 92
			SxxDs[j,k] = cons*Ds[k]*(2*(cb*cb)*FF4 + s2b*FF5 + YB*(c2b*FF6-s2b*FF7));
			SxxDn[j,k] = cons*Dn[k]*(-FF5 + YB*(s2b*FF6 + c2b*FF7));
			SyyDs[j,k] = cons*Ds[k]*(2*(sb*sb)*FF4 - s2b*FF5 - YB*(c2b*FF6-s2b*FF7));
			SyyDn[j,k] = cons*Dn[k]*(-FF5 - YB*(s2b*FF6 + c2b*FF7));
			SxyDs[j,k] = cons*Ds[k]*(s2b*FF4 - c2b*FF5 + YB*(s2b*FF6+c2b*FF7));
			SxyDn[j,k] = cons*Dn[k]*(-YB*(c2b*FF6 - s2b*FF7));
			
			if HSflag==1
				
				#Add some constants
				A1=(y[j]*c4b - YB*c3b);
				A2=(y[j]*s4b - YB*s3b);

				#  Calculate IMAGE AND SUPPLEMENTAL STRESS components due to unit SHEAR and
				#  NORMAL displacement discontinuity
				SxxiDs = cons*Ds[k]*(FF4i - 3*(c2b*FF4i - s2b*FF5i) +
				(2*y[j]*(cb - 3*c3b) + 3*YB*c2b)*FF6i +
				(2*y[j]*(sb - 3*s3b) + 3*YB*s2b)*FF7i -
				2*y[j]*A1*FF8i -
				2*y[j]*A2*FF9i);

				SxxiDn = cons*Dn[k]*(FF5i + (2*y[j]*(sb - 2*s3b) +
				3*YB*s2b)*FF6i - (2*y[j]*(cb - 2*c3b) +
				3*YB*c2b)*FF7i - 2*y[j]*A2*FF8i +
				2*y[j]*A1*FF9i);

				SyyiDs = cons*Ds[k]*(FF4i - (c2b*FF4i - s2b*FF5i) -
				(4*y[j]*sb*s2b - YB*c2b)*FF6i +
				(4*y[j]*sb*c2b + YB*s2b)*FF7i +
				2*y[j]*A1*FF8i +
				2*y[j]*A2*FF9i);

				SyyiDn = cons*Dn[k]*(FF5i - (2*y[j]*sb - YB*s2b)*FF6i +
				(2*y[j]*cb - YB*c2b)*FF7i +
				2*y[j]*A2*FF8i -
				2*y[j]*A1*FF9i);

				SxyiDs = cons*Ds[k]*(s2b*FF4i + c2b*FF5i +
				(2*y[j]*sb*(1+4*c2b) - YB*s2b)*FF6i +
				(2*y[j]*cb*(3-4*c2b) + YB*c2b)*FF7i +
				2*y[j]*A2*FF8i -
				2*y[j]*A1*FF9i);

				SxyiDn = cons*Dn[k]*((4*y[j]*sb*s2b + YB*c2b)*FF6i -
				(4*y[j]*sb*c2b - YB*s2b)*FF7i -
				2*y[j]*A1*FF8i -
				2*y[j]*A2*FF9i);
				
				SxxDs[j,k]=SxxDs[j,k]+SxxiDs;
				SxxDn[j,k]=SxxDn[j,k]+SxxiDn;
				SyyDs[j,k]=SyyDs[j,k]+SyyiDs;
				SyyDn[j,k]=SyyDn[j,k]+SyyiDn;
				SxyDs[j,k]=SxyDs[j,k]+SxyiDs;
				SxyDn[j,k]=SxyDn[j,k]+SxyiDn;
				
			end #half space
			
		end #stress components
		
		if DispFlag==1
			#Calculate the displacement components using eqs. 5.5.4 of C&S, p. 91.
			UxDs[j,k] = Ds[k]*(-pr1*sb*FF2 + pr2*cb*FF3 + YB.*(sb*FF4 - cb*FF5));
			UxDn[j,k] = Dn[k]*(-pr1*cb*FF2 - pr2*sb*FF3 - YB.*(cb*FF4 + sb*FF5));
			UyDs[j,k] = Ds[k]*(+pr1*cb*FF2 + pr2*sb*FF3 - YB.*(cb*FF4 + sb*FF5));
			UyDn[j,k] = Dn[k]*(-pr1*sb*FF2 + pr2*cb*FF3 - YB.*(sb*FF4 - cb*FF5));	
			
			if HSflag==1;
			
				#Add some constants
				A1=(y[j]*c3b - YB*c2b);
				A2=(y[j]*s3b - YB*s2b);
				A3=(y[j]*s2b - YB*sb);
				A4=(y[j]*c2b - YB*cb);
			
				#  See equations 7.4.8 and 7.4.9 in Crouch and Starfield
				#  Calculate image and supplemental displacement components due to unit shear displacement discontinuity
				Uxi_s =Ds[k]* (pr1*sb*FF2i - 
				pr2*cb*FF3i +
				(pr3*A3 + 2*y[j]*s2b)*FF4i +
				(pr3*A4 - y[j]*(1-2*c2b))*FF5i +
				2*y[j]*A2*FF6i -
				2*y[j]*A1*FF7i);

				Uyi_s =Ds[k]* (-pr1*cb*FF2i - 
				pr2*sb*FF3i-
				(pr3*A4 + y[j]*(1-2*c2b))*FF4i+  
				(pr3*A3 - 2*y[j]*s2b)*FF5i+	 
				2*y[j]*A1*FF6i+
				2*y[j]*A2*FF7i);     

				#Calculate image and supplemental displacement components due to unit normal displacement discontinuity
				Uxi_n =Dn[k]*(pr1*cb*FF2i + 
				pr2*sb*FF3i-
				(pr3*A4 - y[j])*FF4i+
				pr3*A3*FF5i-
				2*y[j]*A1*FF6i-
				2*y[j]*A2*FF7i);    

				Uyi_n =Dn[k]* (pr1*sb*FF2i - 
				pr2*cb*FF3i-
				pr3*A3*FF4i-
				(pr3*A4 + y[j])*FF5i+
				2*y[j]*A2*FF6i-
				2*y[j]*A1*FF7i);
				
				UxDs[j,k]=UxDs[j,k]+Uxi_s;
				UxDn[j,k]=UxDn[j,k]+Uxi_n;
				UyDs[j,k]=UyDs[j,k]+Uyi_s;	
				UyDn[j,k]=UyDn[j,k]+Uyi_n;
				
			end #halfspace
			
		end #disp components
		
	end #for every el 
	
end #for every point in space
	
if DispFlag==0
	return(SxxDs,SyyDs,SxyDs,SxxDn,SyyDn,SxyDn) #sum outside when needed
	
elseif StressFlag==0
	return(UxDs,UxDn,UyDs,UyDn)
	
else
	return(SxxDs,SyyDs,SxyDs,SxxDn,SyyDn,SxyDn,UxDs,UxDn,UyDs,UyDn) #sum outside when needed
	
end

end #End LD


function AbbreviateTerms(YB,XB,a)
# First abbreviate repeated terms in the derivatives of f(x,y):
Y2 = YB^2;
XMa = XB-a; XPa = XB+a;
XMa2 = XMa^2; XPa2 = XPa^2;
R1S = XMa2 + Y2; R1S2 = R1S^2;
R2S = XPa2 + Y2; R2S2 = R2S^2;
return(Y2,XMa,XPa,XMa2,XPa2,R1S,R2S,R1S2,R2S2)
end

function CompFF2AndFF3(R1S,R2S,con,XB,YB,XMa,XPa,a)
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
return(FF2,FF3)
end

function CompFF4ToFF7(R1S,R2S,YB,XMa,XPa,XMa2,XPa2,R1S2,R2S2,Y2,con)
FF4 = con*(YB/R1S - YB/R2S);
FF5 = con*(XMa/R1S - XPa/R2S);
# The following derivatives are eqs 553a and b of C&S, p 91
FF6 = con*((XMa2 - Y2)/R1S2 - (XPa2 - Y2)/R2S2);
FF7 = 2*con*YB*(XMa/R1S2 - XPa/R2S2);
return(FF4,FF5,FF6,FF7)
end