# usage #1:
# [Xnw,Ynw] = RotateObject2D!(X,Y,Pa,Pb,Ax1,Ax2)
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
# Ct            - Rotation (defined away from X axis). Direction cosine (CosAx)
#
# St            - (SinAx). Direction cosine (CosAy)
#
# Pa            - Current centre of object at 'X'. Used to centre object at origin before rotation. Use single number
#
# Pb            - Current centre of object at 'Y'. Used to centre object at origin before rotation. Use single number
#
# Arguments: (output)
# X,Y   - the new X,Y point values (overwrites old vals) 
#
#  Author: Tim Davis
#  Copyright 2018, Tim Davis, Potsdam University

@inline @fastmath function RotateObject2D!(X::Float64,Y::Float64,Pa,Pb,Ct,St)
#For single values
#Rotate to new axes Ax Ay Az
x =(Ct[1]*(X-Pa))+(-St[1]*(Y-Pb));
Y =(St[1]*(X-Pa))+( Ct[1]*(Y-Pb));
X=x;
#y=Y;
return(X,Y)
end

@inline @fastmath function RotateObject2D!(X::Array,Y::Array,Pa,Pb,Ct,St)
@simd for i=1:length(X)
	#Rotate to new axes Ax Ay Az
	x=   (Ct[1]*(X[i]-Pa))+(-St[1]*(Y[i]-Pb));
	Y[i]=(St[1]*(X[i]-Pa))+( Ct[1]*(Y[i]-Pb));
	X[i]=x;
	#Y[i]=y;
end

#Vectorised form of: Eq 2.23, Pollard, arranging cosines of new directions in table
#http://continuummechanics.org/stressxforms.html
return(X,Y)
end

