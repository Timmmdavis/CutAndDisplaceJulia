function RotateObject2D(X,Y,Pa,Pb,Ct,St)
# usage #1:
# [Xnw,Ynw] = RotateObject3dNewCoords(X,Y,Pa,Pb,Ax1,Ax2)
#
# The transpose of
#
# A=[Ct,-St;
#    St,Ct];
#
# (i.e.., A') will transform the coordinates from X1X2 into x1x2.
#
# Arguments: (input)
# X             - list of x points (mat or vect)
#
# Y             - list of y points (mat or vect)
#
# Ax1           - The new orientation of the x-axis  (direction cosine 1*3)
#
# Ax2           - The new orientation of the y-axis  (direction cosine 1*3)
#
# Pa            - Centre object 'XY' relative to point PaPbPc
#
# Pb            - Centre object 'XY' relative to point PaPbPc
#
# Arguments: (output)
# x,y   - the new X,Y point values. 
#
#  Author: Tim Davis
#  Copyright 2018, Tim Davis, Potsdam University


x = Array{Float64}(undef, length(X),1);
y = Array{Float64}(undef, length(X),1);
Col = Array{Float64}(undef, 2,1);

for i=1:length(X) #For every point in space

	#Centre object 'XYZ' relative to point PaPbPc
	Col[1] =  X[i]-Pa;
	Col[2] =  Y[i]-Pb;
	
	#Rotate to new axes Ax Ay Az
	x[i]=(Ct*Col[1])+(-St*Col[2]);
	y[i]=(St*Col[1])+( Ct*Col[2]);
	
	#Vectorised form of: Eq 2.23, Pollard, arranging cosines of new directions in table
	#http://continuummechanics.org/stressxforms.html
		
end

return(x,y)

end

