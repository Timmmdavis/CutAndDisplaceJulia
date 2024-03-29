#Test case comparing Okada dislocations with TDs (displacement)

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
(StrikeSlipCosine,DipSlipCosine,FaceNormalVector)=CutAndDisplaceJulia.CreateDirectionCosinesFromStrikeAndDip(Strike,Dip)
#Now rotate the flat plane to the correct location. 
(X,Y,Z)=CutAndDisplaceJulia.RotateObject3DNewCoords(X,Y,Z,0,0,0,DipSlipCosine,StrikeSlipCosine,FaceNormalVector)
Drop=maximum(Z);
Z=Z.-Drop.-TipDepth; #Move fault down
P1=[X[1] Y[1] Z[1];X[2] Y[2] Z[2]] #const 
P2=[X[5] Y[5] Z[5];X[4] Y[4] Z[4]] #const 
P3=[X[3] Y[3] Z[3];X[6] Y[6] Z[6]] #const 

#Repeating array if you want to test with more tris (not the solution vs okada, just speed)
#P1=repeat(P1,50,1)
#P2=repeat(P2,50,1)
#P3=repeat(P3,50,1)
DssVec=repeat([Dss],size(P1,1),1) #const 
DdsVec=repeat([Dds],size(P1,1),1) #const 
DnVec=repeat([Dn],size(P1,1),1) #const 

#= draw fault
using PyPlot
scatter(X,Y,abs.(Z),Z)
cbar = colorbar()
=#

# Start some vectors (spaced points)
xx = range(-10,stop=10,length=50); #linspace deprecated
(xx,yy)=CutAndDisplaceJulia.meshgrid(xx,xx);
zz=ones(size(xx))*-0; #Ground surface

#Get lengths (for reshapes later)
dimx,dimy = size(xx);
#Turn to col vectors
x=reshape(xx,length(xx),1); #const 
y=reshape(yy,length(yy),1); #const 
z=reshape(zz,length(zz),1); #const 

DispFlag=1; #const 
StrainFlag=1; #const 
HSflag=1; #const 

G=ShearModulus(1.); 
ν=PoissonsRatio(0.25);
(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(G,ν);

TotalSlip=sqrt(Dds.^2+Dss.^2);
Rake=90-atand(Dss/Dds); #degrees

println("Vars created -> to TD func")

#using BenchmarkTools
#@btime (No output when you use it)

StrainInfVector=Strains([],[],[],[],[],[]);
DispInfVector=Disps([],[],[]);
@time (StrainAtInfPoints,DispAtInfPoints)= 
CutAndDisplaceJulia.TD(x,y,z,P1,P2,P3,DssVec,DdsVec,DnVec,ν,G,DispFlag,StrainFlag,HSflag,StrainInfVector,DispInfVector)
 
Ux=DispAtInfPoints.Ux;
Uy=DispAtInfPoints.Uy;
Uz=DispAtInfPoints.Uz;

Exx=StrainAtInfPoints.εxx;
Eyy=StrainAtInfPoints.εyy;
Ezz=StrainAtInfPoints.εzz;
Exy=StrainAtInfPoints.εxy;
Exz=StrainAtInfPoints.εxz;
Eyz=StrainAtInfPoints.εyz;  
#= (Exx,Eyy,Ezz,Exy,Exz,Eyz,Ux,Uy,Uz)=
 CutAndDisplaceJulia.TD_sum(ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
 ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
	   ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,
	   UxDn,UyDn,UzDn,
	   UxDss,UyDss,UzDss,
       UxDds,UyDds,UzDds)=#

println("Out of func, too Okada")

MidDeptha=sind(Dip)*Width;
MidDepth=MidDeptha/2+TipDepth;
println("Into Okada func")
(uX,uY,uZ,exx,eyy,ezz,exy,exz,eyz)=CutAndDisplaceJulia.Okada1985RectangularDislocation(x,y,MidDepth,Strike,Dip,Length,Width,Rake,TotalSlip,Dn,ν);

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






