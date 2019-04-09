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
	for i=1:length(Y)
	    X[i]=ObjFunc(Y[i]);
	end

	#@info X Y Desired_X
	
	##Remove 0 before interpolation
	Y=convert(Array,Y)#So we can work on it
	Indx=findall(X.==0.0);
	X=X[length(Indx):end]
	Y=Y[length(Indx):end]

	#using Interpolations
	itp = interpolate((X,),Y, Gridded(Linear()))
	Yv=itp(Desired_X)
	#MATLAB way:
	#Y_Got = interp1(X,Y,X_Desired,'linear');
	@info Yv

	return Yv
end