#Test case comparing Okada dislocations with LDs (displacement). Checks the half space displacement eqs are good. 

#Start creating vars for function: 
println("creating func vars")

#Inputs

#Fault geom
Strike=0;
Dip=45;
Length=100000;
Width=5; # fault width in the DIP direction (WIDTHTD > 0)
TipDepth=1;
Rake=90;
#Fault slip
Dds=0.1;
Dn=0.1;
Dss=0;
#Elastic cons
G=ShearModulus(5.); 
ν=PoissonsRatio(0.25);
(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(G,ν);

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

(SxxDs,SyyDs,SxyDs,SxxDn,SyyDn,SxyDn,UxDs,UxDn,UyDs,UyDn)=CutAndDisplaceJulia.LD(x,z,0,-MidDepth,a,Beta,Ds,Dn,ν,G,DispFlag,StressFlag,HSflag)
#Accumulating arrays
(sXX,sZZ,sXZ,uX,uZ)=CutAndDisplaceJulia.LD_sum(SxxDs,SxxDn,SyyDs,SyyDn,SxyDs,SxyDn,UxDs,UxDn,UyDs,UyDn)
(eXX,eZZ,eXZ)=CutAndDisplaceJulia.HookesLaw2dStress2Strain( sXX,sZZ,sXZ,E,ν )
println("Out of func, too Okada")


println("Into Okada func")
(Ux,Uy,Uz,exx,eyy,ezz,exy,exz,eyz)=CutAndDisplaceJulia.Okada1985RectangularDislocation(x,y,MidDepth,Strike,Dip,Length,Width,Rake,Dds,Dn,ν);
(sxx,szz,sxz)=CutAndDisplaceJulia.HookesLaw2dStrain2Stress( exx,ezz,exz,E,ν,G )
println("Out of func, drawing time, start by reshape")

UxRes=maximum(Ux[:].-uX[:]);
UzRes=maximum(Uz[:].-uZ[:]);

ExxRes=maximum(exx[:].-eXX[:]);
EzzRes=maximum(ezz[:].-eZZ[:]); #Ezz in reality
ExzRes=maximum(exz[:].-eXZ[:]); #Ezz in reality


#using PyPlot
#close()
#scatter(x,sXZ)
#scatter(x,sxz)

#Check LD solution is also traction free
maxSxx=maximum(sXX);
maxSzz=maximum(sZZ);
maxSxz=maximum(sXZ);

println("Values of residuals: LD-Halfspace vs Okada")
@info UxRes UzRes ExxRes EzzRes ExzRes maxSxx maxSzz maxSxz  #Display values in test output
if UxRes>1E-7
	error("UxRes too high, Okada and LD not matching for ground displacement")
end
if UzRes>1E-7
	error("UzRes too high, Okada and LD not matching for ground displacement")
end
if ExxRes>1E-7
	error("ExxRes too high, Okada and LD not matching for ground displacement")
end
if EzzRes>1E-7
	error("EzzRes too high, Okada and LD not matching for ground displacement")
end
if ExzRes>1E-7
	error("ExzRes too high, Okada and LD not matching for ground displacement")
end
if maxSzz>1E-16 || maxSxz>1E-16
	error("LD solution has tractions on the free surface")
end

println("Test Passed")

