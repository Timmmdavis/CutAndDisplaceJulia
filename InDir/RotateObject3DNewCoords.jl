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
Col = Array{Float64}(undef, 3,1);

for i=1:length(X) #For every point in space

	#Centre object 'XYZ' relative to point PaPbPc
	Col[1] =  X[i]-Pa;
	Col[2] =  Y[i]-Pb;
	Col[3] =  Z[i]-Pc;

	#Rotate to new axes Ax Ay Az
	x[i]=(Ax1[1]*Col[1])+(Ax2[1]*Col[2])+(Ax3[1]*Col[3]);
	y[i]=(Ax1[2]*Col[1])+(Ax2[2]*Col[2])+(Ax3[2]*Col[3]);
	z[i]=(Ax1[3]*Col[1])+(Ax2[3]*Col[2])+(Ax3[3]*Col[3]);

	
end


return(x,y,z)

end

