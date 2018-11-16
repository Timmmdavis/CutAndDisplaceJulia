function RotateObject3DNewCoords_prealloc!(ret::AbstractArray{x},X,Y,Z,Pa,Pb,Pc,Ax1,Ax2,Ax3) where x
# usage #1:
#x = Array{Float64}(undef, length(X),1);
#y = Array{Float64}(undef, length(X),1);
#z = Array{Float64}(undef, length(X),1);
#RotateObject3DNewCoords_prealloc!(x,X,Y,Z,P2_1,P2_2,P2_3,Vx[1],Vy[1],Vz[1])
#RotateObject3DNewCoords_prealloc!(y,X,Y,Z,P2_1,P2_2,P2_3,Vx[2],Vy[2],Vz[2])
#RotateObject3DNewCoords_prealloc!(z,X,Y,Z,P2_1,P2_2,P2_3,Vx[3],Vy[3],Vz[3])
#  Author: Tim Davis
#  Copyright 2018, Tim Davis, Potsdam University

for i in 1:length(X) #For every point in space

	#Rotate to new axes Ax Ay Az
	ret[i]=(Ax1[1]*(X[i]-Pa))+(Ax2[1]*(Y[i]-Pb))+(Ax3[1]*(Z[i]-Pc));
	#For x give first row of x, for y 2nd etc. 
		
end

nothing 
end

