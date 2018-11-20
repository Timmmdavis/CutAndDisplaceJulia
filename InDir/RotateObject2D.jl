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
# Ct            - Rotation (defined away from X axis). Direction cosine (CosAx)
#
# St            - (SinAx). Direction cosine (CosAy)
#
# Pa            - Current centre of object at 'X'. Used to centre object at origin before rotation. Use single number
#
# Pb            - Current centre of object at 'Y'. Used to centre object at origin before rotation. Use single number
#
# Arguments: (output)
# x,y   - the new X,Y point values. 
#
#  Author: Tim Davis
#  Copyright 2018, Tim Davis, Potsdam University

#Only move coords if needed
if Pa!=0 || Pb!=0 
	X =  X.-Pa;
	Y =  Y.-Pb;
end

#Rotate to new axes Ax Ay Az
x=(Ct[1].*X).+(-St[1].*Y);
y=(St[1].*X).+( Ct[1].*Y);

#Vectorised form of: Eq 2.23, Pollard, arranging cosines of new directions in table
#http://continuummechanics.org/stressxforms.html


return(x,y)

end

