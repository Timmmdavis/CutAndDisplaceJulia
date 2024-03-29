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
	lps=0
	max_reached=false
	for i=1:length(Y)
		lps+=1
	    X[i]=ObjFunc(Y[i]);
	    #Stops bad knot vector error (sometimes the loops in the obj func are not enough to stop neg values coming out)
	    if X[i]<0
	    	X[i]=1e-12*i
	    end
	    #break loop early if we have already hit DesiredVal
	    if X[i]>Desired_X
	    	max_reached=true
	    	X=X[1:i];
	    	Y=Y[1:i];
	    	break
	    end
	end

	
	#if we broke out of first loop without reaching max increase max pressure
	if max_reached==false
		println("Didnt reach max - retrying")
		Y = range(MaxVal,stop=MaxVal*2,length=NumberOfIts); #linspace deprecated
		X=zeros(size(Y));
		breakat=[]
		lps=0
		max_reached=false
		for i=1:length(Y)
			lps+=1
		    X[i]=ObjFunc(Y[i]);
		    #Stops bad knot vector error (sometimes the loops in the obj func are not enough to stop neg values coming out)
		    if X[i]<0
		    	X[i]=1e-12*i
		    end
		    #break loop early if we have already hit DesiredVal
		    if X[i]>Desired_X
		    	max_reached==true
		    	X=X[1:i];
		    	Y=Y[1:i];
		    	break
		    end
		end
	end

	#if we broke out of last loop super early redo with a better sampling 
	Minimum=3
	if lps<Minimum 
		println("droppingdownMaxPressure")
		if length(Y)<=2
			Y = range(MinVal,maximum(Y),length=NumberOfIts); #linspace deprecated
		else
			Y = range(MinVal,stop=Y[Minimum],length=NumberOfIts); #linspace deprecated
		end
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
		    	@bp
		    	max_reached==true
		    	X=X[1:i];
		    	Y=Y[1:i];
		    	break
		    end
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