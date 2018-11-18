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
#Only move coords if needed
if Pa!=0 && Pb!=0 
	if length(X)==1 #Catch error if X is not an array
		X =  X-Pa;
		Y =  Y-Pb;
	else
		for i=1:length(X) #And if it is an array do for every point 
			#Centre object 'XYZ' relative to point PaPbPc
			X[i] =  X[i]-Pa;
			Y[i] =  Y[i]-Pb;
		end
	end
end

for i=1:length(X) #For every point in space
	
	#Rotate to new axes Ax Ay Az
	x[i]=(Ct[1]*X[i])+(-St[1]*Y[i]);
	y[i]=(St[1]*X[i])+( Ct[1]*Y[i]);
	
	#Vectorised form of: Eq 2.23, Pollard, arranging cosines of new directions in table
	#http://continuummechanics.org/stressxforms.html
		
end

return(x,y)

end

