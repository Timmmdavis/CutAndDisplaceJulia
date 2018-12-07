#Test case comparing Okada dislocations with TDs (displacement)

#Start creating vars for function: 
println("creating func vars")

#Inputs

#Fault geom
Strike=90;
Dip=45;
Length=1000;
Width=5; # fault width in the DIP direction (WIDTH > 0)
TipDepth=1;
Rake=90;
#Fault slip
Dds=1;
Dn=0;
Dss=0;
#Elastic cons
nu=0.25;

#Dip relative to X-axis
Beta=rad2deg(90-Dip);
a=Width;
#Depth of midpoint
MidDeptha=sind(Dip)*Width;
MidDepth=MidDeptha/2+TipDepth;


# Start some vectors (spaced points)
x = range(-10,stop=10,length=50); #linspace deprecated
(x,y)=MyModule.meshgrid(x,x);
z=ones(size(x))*-0; #Ground surface

#Get lengths (for reshapes later)
dimx,dimy = size(x);
#Turn to col vectors
x=reshape(x,length(x),1);
y=reshape(y,length(y),1);
z=reshape(z,length(z),1);


DispFlag=1;
StressFlag=1;
HSflag=0;

println("Vars created -> to LD func")

(SxxDs,SyyDs,SxyDs,SxxDn,SyyDn,SxyDn,UxDs,UxDn,UyDs,UyDn)=MyModule.LD(x,y,0,MidDepth,a,Beta,Dds,Dn,nu,mu,DispFlag,StressFlag,HSflag)
#Accumulating arrays
(sXX,sYY,sXY,uX,uY)=MyModule.LD_sum(SxxDs,SxxDn,SyyDs,SyyDn,SxyDs,SxyDn,UxDs,UxDn,UyDs,UyDn)

println("Out of func, too Okada")


println("Into Okada func")
(uX,uY,uZ)=MyModule.Okada1985RectangularDislocation(x,y,MidDepth,Strike,Dip,Length,Width,Rake,Dds,Dn,nu);

println("Out of func, drawing time, start by reshape")


UxRes=maximum(Ux[:].-uX[:]);
UyRes=maximum(Uy[:].-uY[:]);
UzRes=maximum(Uz[:].-uZ[:]);

println("Values of residuals: TDE vs Okada")
@info UxRes UyRes UzRes #Display values in test output
if UxRes>1E-14
	error("UxRes too high, Okada and TD not matching for displacement")
end
if UyRes>1E-14
	error("UyRes too high, Okada and TD not matching for displacement")
end
if UzRes>1E-14
	error("UzRes too high, Okada and TD not matching for displacement")
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

