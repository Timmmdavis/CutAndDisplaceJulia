function TestFooAllocs(x,y,z,p1,p2,p3,Ct,St,Ax1,Ax2,Ax3)


# @time (X2,Y2)=MyModule.RotateObject2D(x,y,p1,p2,Ct,St)
# @time (x,y)= MyModule.RotateObject2D!(x,y,p1,p2,Ct,St)


# if any(x.!=X2) || any(y.!=Y2)
	# error("you have introduced an error")
# end
# end


@time (X2,Y2,Z2)= MyModule.RotateObject3DNewCoords(x,y,z,p1,p2,p3,Ax1,Ax2,Ax3)
@time (x,y,z)= MyModule.RotateObject3DNewCoords!(x,y,z,p1,p2,p3,Ax1,Ax2,Ax3)

if any(x.!=X2) || any(y.!=Y2) || any(z.!=Z2)
	error("you have introduced an error")
end
end