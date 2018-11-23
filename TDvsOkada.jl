#Start creating vars for function: 
println("creating func vars")

#Define dipping fault plane (Would be better to add rectangle and transform....
P1=[2.59258869000000	0.509504670000000	-5.69846310000000; -2.59258869000000	-0.509504670000000	-1]
P2=[1.73753833000000	1.99049533000000	-1; -1.73753833000000	-1.99049533000000	-5.69846310000000]
P3=[-1.73753833000000	-1.99049533000000	-5.69846310000000; 1.73753833000000	1.99049533000000	-1]
using PyPlot
plot_trisurf(P1,P2,P3);


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
