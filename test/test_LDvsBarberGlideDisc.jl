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
(X,Y) = CutAndDisplaceJulia.meshgrid(x,x);
#Get lengths (for reshapes later)
dimx,dimy = size(X);

#Doing with two Els
a=1;
hlflen = [0.5 0.5];  
xe=[-0.5 0.5];
ye=[0. 0.];
Beta=[0. 0.];
Ds=[1. 1.];
Dn=[0. 0.];


DispFlag=1;
StressFlag=1;
HSflag=0;

println("Vars created -> to LD func")

(SxxDs,SyyDs,SxyDs,SxxDn,SyyDn,SxyDn,UxDs,UxDn,UyDs,UyDn)=CutAndDisplaceJulia.LD(X[:],Y[:],xe,ye,hlflen,Beta,Ds,Dn,nu,mu,DispFlag,StressFlag,HSflag)
#Accumulating arrays
(sXX,sYY,sXY,uX,uY)=CutAndDisplaceJulia.LD_sum(SxxDs,SxxDn,SyyDs,SyyDn,SxyDs,SxyDn,UxDs,UxDn,UyDs,UyDn)

println("Out of func, too Barber1992_GlideDislocation")

(X,Y,Sxx,Syy,Sxy,Ux,Uy)=CutAndDisplaceJulia.Barber1992_GlideDislocation(k,mu,X[:],Y[:],a,Ds[1],nu)


println("Out of func, drawing time, start by reshape")

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

# ##= if you want to draw remove this line
# x=reshape(X,dimx,dimy);
# y=reshape(Y,dimx,dimy);
# uX=reshape(uX,dimx,dimy);
# Ux=reshape(Ux,dimx,dimy);
# uY=reshape(uY,dimx,dimy);
# Uy=reshape(Uy,dimx,dimy);

# Value=uX
# using NaNMath
# Top=maximum([NaNMath.maximum(Value),abs(NaNMath.minimum(Value))])
# steps=10; #Steps from centre to top. 
# levels = [-Top:Top/steps:Top;]
# using PyPlot
# close()
# contourf(x,y,Value, levels=levels);
# cbar = colorbar()
# ##=#

