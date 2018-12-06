#Start creating vars for function: 
println("creating func vars")

#Inputs
nu = 0.25;
k = 3-(4*nu); 
mu = 500;

#create x and y
spacing=0.1;
minx=-4; maxx=4;
x = [minx:spacing:maxx;];
(X,Y) = MyModule.meshgrid(x,x);

a = 1;  
b=0.0001; 
Beta=0;
Ds=1;
Dn=0.;

println("Vars created -> to LD func")

(uX,uY)=MyModule.LDdispFS(X[:],Y[:],0,0,a,Beta,Ds,Dn,nu)
(sXX,sYY,sXY)=MyModule.LDstressFS(X[:],Y[:],0,0,a,Beta,Ds,Dn,nu,mu)

println("Out of func, too Barber1992_GlideDislocation")

(X,Y,Sxx,Syy,Sxy,Ux,Uy)=MyModule.Barber1992_GlideDislocation(k,mu,X[:],Y[:],a,Ds,nu)


println("Out of func, drawing time, start by reshape")

@info Ux[1:10] uX[1:10]

UxRes=maximum(filter(!isnan, Ux[:].-uX[:]));
UyRes=maximum(filter(!isnan, Uy[:].-uY[:]));
SxxRes=maximum(filter(!isnan, Sxx[:].-sXX[:]));
SyyRes=maximum(filter(!isnan, Syy[:].-sYY[:]));
SxyRes=maximum(filter(!isnan, Sxy[:].-sXY[:]));

println("Values of residuals: LD vs Barber")
@info UxRes UyRes SxxRes SyyRes SxyRes
if UxRes>1E-13
	error("UxRes too high, LD and glide dislocation not matching for displacement")
end
if UyRes>1E-13
	error("UyRes too high, LD and glide dislocation not matching for displacement")
end
if SxxRes>1E-12
	error("SxxRes too high, LD and glide dislocation not matching for stress")
end
if SyyRes>1E-12
	error("SyyRes too high, LD and glide dislocation not matching for stress")
end
if SxyRes>1E-12
	error("SxyRes too high, LD and glide dislocation not matching for stress")
end


println("Test Passed")

# #= if you want to draw remove this line
# x=reshape(x,dimx,dimy);
# y=reshape(y,dimx,dimy);
# uX=reshape(uX,dimx,dimy);
# Ux=reshape(Ux,dimx,dimy);
# uY=reshape(uY,dimx,dimy);
# Uy=reshape(Uy,dimx,dimy);
# uZ=reshape(uZ,dimx,dimy);
# Uz=reshape(Uz,dimx,dimy);


# Value=uZ
# using NaNMath
# Top=maximum([NaNMath.maximum(Value),abs(NaNMath.minimum(Value))])
# steps=10; #Steps from centre to top. 
# levels = [-Top:Top/steps:Top;]
# using PyPlot
# close()
# contourf(x,y,Value, levels=levels);
# cbar = colorbar()
# =#

