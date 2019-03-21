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
StressFlag=1; #const 
HSflag=1; #const 


TotalSlip=sqrt(Dds.^2+Dss.^2);
Rake=90-atand(Dss/Dds); #degrees

G=ShearModulus(1.); 
ν=PoissonsRatio(0.25);
(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(G,ν);

println("Vars created -> to TD func")


#using BenchmarkTools
#@btime (No output when you use it)

@time (ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn,
 ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
 ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,
 UxDn,UyDn,UzDn,
 UxDss,UyDss,UzDss,
 UxDds,UyDds,UzDds)=
 CutAndDisplaceJulia.TD(x,y,z,P1,P2,P3,DssVec,DdsVec,DnVec,ν,G,DispFlag,StressFlag,HSflag)
 
#Converting this to stress tensor influences. 
#Uses Hooke's Law to convert strain to stress. Equation 7.131 and 7.132 in David Pollards Book. 
(SxxDss,SyyDss,SzzDss,SxyDss,SxzDss,SyzDss) = CutAndDisplaceJulia.HookesLaw3dStrain2Stress(ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,λ,G);
(SxxDds,SyyDds,SzzDds,SxyDds,SxzDds,SyzDds) = CutAndDisplaceJulia.HookesLaw3dStrain2Stress(ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,λ,G);
(SxxDn,SyyDn,SzzDn,SxyDn,SxzDn,SyzDn) = CutAndDisplaceJulia.HookesLaw3dStrain2Stress(ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn,λ,G);


CosAx=FaceNormalVector[:,1];
CosAy=FaceNormalVector[:,2];
CosAz=FaceNormalVector[:,3];
CutAndDisplaceJulia.CreateTriangleNormal(P1,P2,P3)

#error("CheckOutput")
#DssTn=CutAndDisplaceJulia.CalculateNormalTraction3d( SxxDss,SyyDss,SzzDss,SxyDss,SxzDss,SyzDss,CosAx,CosAy,CosAz )
#DdsTn=CutAndDisplaceJulia.CalculateNormalTraction3d( SxxDds,SyyDds,SzzDds,SxyDds,SxzDds,SyzDds,CosAx,CosAy,CosAz )
#DnTn =CutAndDisplaceJulia.CalculateNormalTraction3d( SxxDn,SyyDn,SzzDn,SxyDn,SxzDn,SyzDn,CosAx,CosAy,CosAz )
