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
# Pa            - Current centre of object at 'X'. Used to centre object at origin before rotation. Use single number
#
# Pb            - Current centre of object at 'Y'. Used to centre object at origin before rotation. Use single number
#
# Pc            - Current centre of object at 'Z'. Used to centre object at origin before rotation. Use single number
#
# Arguments: (output)
# x,y,z   - the new X,Y and Z point values. 
#
#
# X=[0 0 1 1 0 0 0 0 0 1 1 0 1 1 1 1];
# Y=[0 0 0 0 0 1 1 0 1 1 1 1 1 0 0 1];
# Z=[0 1 1 0 0 0 1 1 1 1 0 0 0 0 1 1];
# X=convert(Array{Float64,2},X)
# Y=convert(Array{Float64,2},Y)
# Z=convert(Array{Float64,2},Z)
# using PyPlot
# plot_wireframe(X,Y,Z, color="red")
# #flip Y and Z axes
# (x,y,z)=CutAndDisplaceJulia.RotateObject3DNewCoords(X,Y,Z,0,0,0,[1;0;0],[0;-1;0],[0;0;1])
# plot_wireframe(x,y,z, color="blue")
#
#  Author: Tim Davis
#  Copyright 2018, Tim Davis, Potsdam University

#Only move coords if needed
if Pa!=0 || Pb!=0 || Pc!=0
	X =  X.-Pa;
	Y =  Y.-Pb;
	Z =  Z.-Pc;	
end

#Rotate to new axes Ax Ay Az
x=(Ax1[1].*X).+(Ax2[1].*Y).+(Ax3[1].*Z);
y=(Ax1[2].*X).+(Ax2[2].*Y).+(Ax3[2].*Z);
z=(Ax1[3].*X).+(Ax2[3].*Y).+(Ax3[3].*Z);	

return(x,y,z)

end
