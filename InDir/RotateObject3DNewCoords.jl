function RotateObject3DNewCoords(X,Y,Z,Pa,Pb,Pc,Ax1,Ax2,Ax3)
# usage #1:
# [Xnw,Ynw,Znw] = RotateObject3dNewCoords(X,Y,Z,Pa,Pb,Pc,Ax1,Ax2,Ax3)
#
# The transpose of A=[Ax;Ay;Az] (i.e.., A') will transform the coordinates from X1X2X3 into x1x2x3.
#
# Arguments: (input)
# X             - list of x points (mat or vect)
#
# Y             - list of y points (mat or vect)
#
# Z             - list of y points (mat or vect)
#
# Ax1           - The new orientation of the x-axis  (direction cosine 1*3)
#
# Ax2           - The new orientation of the y-axis  (direction cosine 1*3)
#
# Ax3           - The new orientation of the z-axis  (direction cosine 1*3)
#
# Pa            - Centre object 'XYZ' relative to point PaPbPc
#
# Pb            - Centre object 'XYZ' relative to point PaPbPc
#
# Pc            - Centre object 'XYZ' relative to point PaPbPc
#
# Arguments: (output)
# x,y,z   - the new X,Y and Z point values. 
#
#  Author: Tim Davis
#  Copyright 2018, Tim Davis, Potsdam University

x = Array{Float64}(undef, length(X),1);
y = Array{Float64}(undef, length(X),1);
z = Array{Float64}(undef, length(X),1);

#Only move coords if needed
if Pa!=0 && Pb!=0 && Pc!=0
	if length(X)==1 #Catch error if X is not an array
		X =  X-Pa;
		Y =  Y-Pb;
		Z =  Z-Pc;	
	else
		for i=1:length(X) #And if it is an array do for every point 
			#Centre object 'XYZ' relative to point PaPbPc
			X[i] =  X[i]-Pa;
			Y[i] =  Y[i]-Pb;
			Z[i] =  Z[i]-Pc;
		end
	end
end

for i=1:length(X) #For every point 
	#Rotate to new axes Ax Ay Az
	x[i]=(Ax1[1]*X[i])+(Ax2[1]*Y[i])+(Ax3[1]*Z[i]);
	y[i]=(Ax1[2]*X[i])+(Ax2[2]*Y[i])+(Ax3[2]*Z[i]);
	z[i]=(Ax1[3]*X[i])+(Ax2[3]*Y[i])+(Ax3[3]*Z[i]);	
end


return(x,y,z)

end

