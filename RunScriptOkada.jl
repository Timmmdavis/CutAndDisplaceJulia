#Start creating vars for function: 
println("creating func vars")

# Comment - start some vectors (spaced points)
x = range(-10,stop=10,length=50); #linspace deprecated
(x,y)=MyModule.meshgrid(x,x);
z=ones(size(x))*2;
#Get lengths (for reshapes later)
dimx,dimy = size(x);
#Turn to col vectors
x=reshape(x,length(x),1);
y=reshape(y,length(y),1);
z=reshape(z,length(z),1);


Z=2;
Strike=30;
Dip=70;
Length=5;
Width=3;
Rake=-45;
Ds=1;
Dn=1;
nu=0.25;

println("Vars created -> to func")
(uX,uY,uZ)=MyModule.Okada1985RectangularDislocation(x,y,Z,Strike,Dip,Length,Width,Rake,Ds,Dn,nu);
println("Out of func, drawing time")


Value=uY;

####Some reshaping for drawing:
#Was one func that reshaped stuff
x=reshape(x,dimx,dimy);
y=reshape(y,dimx,dimy);
Value=reshape(Value,dimx,dimy);

#Draw
using NaNMath
Top=maximum([NaNMath.maximum(Value),abs(NaNMath.minimum(Value))])
steps=10; #Steps from centre to top. 
levels = [-Top:Top/steps:Top;]
using PyPlot
close()
contourf(x,y,Value, levels=levels);
cbar = colorbar()