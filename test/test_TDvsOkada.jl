#Test case comparing Okada dislocations with TDs (displacement)


println("ToDo 
		1st Remove coord transform of global sum (Ux etc, do inside loop)
		2nd If this doesnt solve the halfspace FS func issue (allocating as zeros and adding) 
		,then put this before the original call
		3rd clean up funcs, write what each does actually, i.e. coord transform ADCS to TDCS, Comp Vect etc
		4th try and speed up reducing alloc etc")

#Start creating vars for function: 
println("creating func vars")

#Inputs

#Fault geom
Strike=60;
Dip=45;
Length=5;
Width=5;
TipDepth=1;
Rake=90;
#Fault slip
Dds=2.;  #2
Dn=0.6; #0.6
Dss=0.5;#0.5
#Elastic cons
nu=0.25;


###Section - arranging the fault surface for the TDE solution
#Arranged as so
#P1=[-1  1 0 ; 1 -1 0]; #top left and bottom right
#P2=[-1 -1 0 ;-1 -1 0]; #bottom left
#P3=[ 1  1 0 ; 1  1 0]; #top right
X=[-1  1 -1 -1 1 1];
Y=[ 1 -1 -1 -1 1 1];
Z=[ 0  0  0  0 0 0];
X=(X./2).*Width;
Y=(Y./2).*Length;
#Getting the direction cosines for the actual plane
(StrikeSlipCosine,DipSlipCosine,FaceNormalVector)=MyModule.CreateDirectionCosinesFromStrikeAndDip(Strike,Dip)
#Now rotate the flat plane to the correct location. 
(X,Y,Z)=MyModule.RotateObject3DNewCoords(X,Y,Z,0,0,0,DipSlipCosine,StrikeSlipCosine,FaceNormalVector)
Drop=maximum(Z);
Z=Z.-Drop.-TipDepth; #Move fault down
P1=[X[1] Y[1] Z[1];X[2] Y[2] Z[2]]
P2=[X[5] Y[5] Z[5];X[4] Y[4] Z[4]]
P3=[X[3] Y[3] Z[3];X[6] Y[6] Z[6]]

#= draw fault
using PyPlot
scatter(X,Y,abs.(Z),Z)
cbar = colorbar()
=#

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
HSflag=1;
mu=1.;

TotalSlip=sqrt(Dds.^2+Dss.^2);
Rake=90-atand(Dss/Dds); #degrees

#@info TotalSlip Rake
#poop

println("Vars created -> to TD func1")
##3

#using BenchmarkTools
#@btime (No output when you use it)

@time (ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn,
 ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
 ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,
 UxDn,UyDn,UzDn,
 UxDss,UyDss,UzDss,
 UxDds,UyDds,UzDds)=
 MyModule.TD(x,y,z,P1,P2,P3,[Dss Dss],[Dds Dds],[Dn Dn],nu,mu,DispFlag,StressFlag,HSflag)
 
 
 (Exx,Eyy,Ezz,Exy,Exz,Eyz,Ux,Uy,Uz)=
 MyModule.TD_sum(ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
 ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
	   ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,
	   UxDn,UyDn,UzDn,
	   UxDss,UyDss,UzDss,
       UxDds,UyDds,UzDds)

println("Out of func, too Okada")

MidDeptha=sind(Dip)*Width;
MidDepth=MidDeptha/2+TipDepth;
println("Into Okada func")
(uX,uY,uZ,exx,eyy,ezz,exy,exz,eyz)=MyModule.Okada1985RectangularDislocation(x,y,MidDepth,Strike,Dip,Length,Width,Rake,TotalSlip,Dn,nu);

println("Out of func, drawing time, start by reshape")


UxRes=maximum(Ux[:].-uX[:]);
UyRes=maximum(Uy[:].-uY[:]);
UzRes=maximum(Uz[:].-uZ[:]);

ExxRes=maximum(exx[:].-Exx[:]);
EyyRes=maximum(eyy[:].-Eyy[:]);
EzzRes=maximum(ezz[:].-Ezz[:]); 
ExyRes=maximum(exy[:].-Exy[:]);
ExzRes=maximum(exz[:].-Exz[:]);
EyzRes=maximum(eyz[:].-Eyz[:]);

@info UxRes UyRes UzRes ExxRes EyyRes EzzRes ExyRes ExzRes EyzRes  #Display values in test output


println("Values of residuals: TDE vs Okada")
@info UxRes UyRes UzRes #Display values in test output
if UxRes>1E-13
	error("UxRes too high, Okada and TD not matching for displacement")
end
if UyRes>1E-13
	error("UyRes too high, Okada and TD not matching for displacement")
end
if UzRes>1E-13
	error("UzRes too high, Okada and TD not matching for displacement")
end

if ExxRes>1E-13
	error("ExxRes too high, Okada and TD not matching for displacement")
end
if EyyRes>1E-13
	error("EyyRes too high, Okada and TD not matching for displacement")
end
if EzzRes>1E-13
	error("EzzRes too high, Okada and TD not matching for displacement")
end
if ExyRes>1E-13
	error("ExyRes too high, Okada and TD not matching for displacement")
end
if ExzRes>1E-13
	error("ExzRes too high, Okada and TD not matching for displacement")
end
if EyzRes>1E-13
	error("EyzRes too high, Okada and TD not matching for displacement")
end
println("Test P1 Passed")

println("Test P2 Off for now...")
println("Now checking that when we are for only one component e.g. dipslip we only get values in those matricies")
##Checking that there are no runaway components
(ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn,
 ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
 ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,
 UxDn,UyDn,UzDn,
 UxDss,UyDss,UzDss,
 UxDds,UyDds,UzDds)=
 MyModule.TD(x,y,z,P1,P2,P3,[1. 1.],[0. 0.],[0. 0.],nu,mu,DispFlag,StressFlag,HSflag)
uux=UxDn+UxDds;#
uuy=UyDn+UyDds;
uuz=UzDn+UzDds;
eexx=ExxDn+ExxDds;#
eeyy=EyyDn+EyyDds;
eezz=EzzDn+EzzDds;
eexy=ExyDn+ExyDds;
eexz=ExzDn+ExzDds;
eeyz=EyzDn+EyzDds;
Exx=sum(eexx,dims=2);
Eyy=sum(eeyy,dims=2);
Ezz=sum(eezz,dims=2);
Exy=sum(eexy,dims=2);
Exz=sum(eexz,dims=2);
Eyz=sum(eeyz,dims=2);
Ux=sum(uux,dims=2);
Uy=sum(uuy,dims=2);
Uz=sum(uuz,dims=2);
if any((Exx+Eyy+Ezz+Exy+Exz+Eyz+Ux+Uy+Uz).>eps())
	error("Only SS components wanted but its leaking into others.")
end

##Checking that there are no runaway components
(ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn,
 ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
 ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,
 UxDn,UyDn,UzDn,
 UxDss,UyDss,UzDss,
 UxDds,UyDds,UzDds)=
 MyModule.TD(x,y,z,P1,P2,P3,[0. 0.],[1. 1.],[0. 0.],nu,mu,DispFlag,StressFlag,HSflag)
uux=UxDn+UxDss;#
uuy=UyDn+UyDss;
uuz=UzDn+UzDss;
eexx=ExxDn+ExxDss;#
eeyy=EyyDn+EyyDss;
eezz=EzzDn+EzzDss;
eexy=ExyDn+ExyDss;
eexz=ExzDn+ExzDss;
eeyz=EyzDn+EyzDss;
Exx=sum(eexx,dims=2);
Eyy=sum(eeyy,dims=2);
Ezz=sum(eezz,dims=2);
Exy=sum(eexy,dims=2);
Exz=sum(eexz,dims=2);
Eyz=sum(eeyz,dims=2);
Ux=sum(uux,dims=2);
Uy=sum(uuy,dims=2);
Uz=sum(uuz,dims=2);
if any((Exx+Eyy+Ezz+Exy+Exz+Eyz+Ux+Uy+Uz).>eps())
	error("Only DS components wanted but its leaking into others.")
end
 
 ##Checking that there are no runaway components
(ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn,
 ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
 ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,
 UxDn,UyDn,UzDn,
 UxDss,UyDss,UzDss,
 UxDds,UyDds,UzDds)=
 MyModule.TD(x,y,z,P1,P2,P3,[0. 0.],[0. 0.],[1. 1.],nu,mu,DispFlag,StressFlag,HSflag)
uux=UxDss+UxDds;#
uuy=UyDss+UyDds;
uuz=UzDss+UzDds;
eexx=ExxDss+ExxDds;#
eeyy=EyyDss+EyyDds;
eezz=EzzDss+EzzDds;
eexy=ExyDss+ExyDds;
eexz=ExzDss+ExzDds;
eeyz=EyzDss+EyzDds;
Exx=sum(eexx,dims=2);
Eyy=sum(eeyy,dims=2);
Ezz=sum(eezz,dims=2);
Exy=sum(eexy,dims=2);
Exz=sum(eexz,dims=2);
Eyz=sum(eeyz,dims=2);
Ux=sum(uux,dims=2);
Uy=sum(uuy,dims=2);
Uz=sum(uuz,dims=2);
if any((Exx+Eyy+Ezz+Exy+Exz+Eyz+Ux+Uy+Uz).>eps())
	error("Only SS components wanted but its leaking into others.")
end

println("Test P2 Passed")


# ###= if you want to draw remove this line
# x=reshape(x,dimx,dimy);
# y=reshape(y,dimx,dimy);
# Exx=reshape(Exx,dimx,dimy);
# uX=reshape(uX,dimx,dimy);
# Ux=reshape(Ux,dimx,dimy);
# uY=reshape(uY,dimx,dimy);
# Uy=reshape(Uy,dimx,dimy);
# uZ=reshape(uZ,dimx,dimy);
# Uz=reshape(Uz,dimx,dimy);

# Value=Exx
# using NaNMath
# Top=maximum([NaNMath.maximum(Value),abs(NaNMath.minimum(Value))])
# steps=10; #Steps from centre to top. 
# levels = [-Top:Top/steps:Top;]
# using PyPlot
# close()
# contourf(x,y,Value, levels=levels);
# cbar = colorbar()
# #=#






