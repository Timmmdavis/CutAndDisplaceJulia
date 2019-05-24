function WalkAndInterp(ObjFunc, MinVal,MaxVal, NumberOfIts,Desired_X)

	#ObjFunc - pointer to the objective func
	#MaxVal  - Maximum value out of func
	#NumberOfIts - Number of interations
	#DesiredVal - Value we want

	#Walks through an objective function from 0 to the given value and computes
	#the result (Y). Then this interpolates the result for a given X to find Y
	#and returns this.
	#Example : Given volume (DesiredVal) in crack for an increasing pressure
	#(0-MaxVal), we know volume will increase with pressure (unequivocal) and
	#we make sure the ObjFunc accepts pressure as an input and spits out a
	#volume.

	#Step=MaxVal/NumberOfIts;
	Y = range(MinVal,stop=MaxVal,length=NumberOfIts); #linspace deprecated
	X=zeros(size(Y));
	breakat=[]
	for i=1:length(Y)
	    X[i]=ObjFunc(Y[i]);
	    #Stops bad knot vector error (sometimes the loops in the obj func are not enough to stop neg values coming out)
	    if X[i]<0
	    	X[i]=1e-12*i
	    end
	    #break loop early if we have already hit DesiredVal
	    if X[i]>Desired_X
	    	X=X[1:i];
	    	Y=Y[1:i];
	    	break
	    end
	end

	
	#@info X Y Desired_X
	
	#scatter(X,Y,markersize=(Desired_X/25))
	@bp

	#Reduce to only values that surrond our desired X (assuming we stopped one value above this)
	X=X[end-1:end]
	Y=Y[end-1:end]

	#using Interpolations
	itp = interpolate((X,),Y, Gridded(Linear()))
	Yv=TryInterp(itp,Desired_X)

	#@info X 
	#@info Y

	return Yv
end

function TryInterp(itp,Desired_X) 
	try 
		Yv=itp(Desired_X); 
	catch
		println("increase MinVal and MaxVal limits, no extrapolation set")
		Yv=NaN;
	end
end