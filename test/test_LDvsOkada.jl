#Test case comparing Okada dislocations with LDs (displacement). Checks the half space displacement eqs are good. 

#Start creating vars for function: 
println("creating func vars")

#Inputs

#Fault geom
Strike=0;
Dip=0;
Length=100000;
Width=5; # fault width in the DIP direction (WIDTH > 0)
TipDepth=1;
Rake=90;
#Fault slip
Dds=0.1;
Dn=0.1;
Dss=0;
#Elastic cons
nu=0.25;
mu=5;

#Dip relative to X-axis
Beta=deg2rad(180-Dip);
a=Width/2;
Ds=Dds;
#Depth of midpoint
MidDeptha=sind(Dip)*Width;
MidDepth=MidDeptha/2+TipDepth;

# Start some vectors (spaced points)
x = [-10:0.1:10;]
y=zeros(length(x))
z=zeros(length(x))

DispFlag=1;
StressFlag=1;
HSflag=1;

println("Vars created -> to LD func")

(SxxDs,SyyDs,SxyDs,SxxDn,SyyDn,SxyDn,UxDs,UxDn,UyDs,UyDn)=MyModule.LD(x,z,0,-MidDepth,a,Beta,Ds,Dn,nu,mu,DispFlag,StressFlag,HSflag)
#Accumulating arrays
(sXX,sYY,sXY,uX,uZ)=MyModule.LD_sum(SxxDs,SxxDn,SyyDs,SyyDn,SxyDs,SxyDn,UxDs,UxDn,UyDs,UyDn)

println("Out of func, too Okada")


println("Into Okada func")
(Ux,Uy,Uz,unn,une,uen,uee)=MyModule.Okada1985RectangularDislocation(x,y,MidDepth,Strike,Dip,Length,Width,Rake,Dds,Dn,nu);

println("Out of func, drawing time, start by reshape")


UxRes=maximum(Ux[:].-uX[:]);
#UyRes=maximum(Uy[:].-uY[:]);
UzRes=maximum(Uz[:].-uZ[:]);

#close()
#using PyPlot
#scatter(x,Uz)
#scatter(x,uZ)

println("Values of residuals: TDE vs Okada")
@info UxRes UzRes #Display values in test output
if UxRes>1E-7
	error("UxRes too high, Okada and LD not matching for ground displacement")
end
if UzRes>1E-7
	error("UzRes too high, Okada and LD not matching for ground displacement")
end


println("Test Passed")

