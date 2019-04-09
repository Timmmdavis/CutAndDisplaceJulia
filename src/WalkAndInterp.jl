function WalkAndInterp(ObjFunc, MinVal,MaxVal, NumberOfIts,DesiredVal)

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

	#Remove 0 before interpolation
	Indx=find(X.==0);
	Y[Indx[1:end-1]]=[];
	X[Indx[1:end-1]]=[];

	#using Interpolations
	start=minimum(X);
	itp = interpolate(Y, BSpline(Linear()))
	X_Good=itp(DesiredVal+(1-start))
	#MATLAB way:
	#X_Good = interp1(X,Y,DesiredVal,'linear');



	return XGood
end