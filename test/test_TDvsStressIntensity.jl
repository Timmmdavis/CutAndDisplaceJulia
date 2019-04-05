#Test case comparing to Penny shaped crack

#Start creating vars for function: 
println("creating func vars")

#Load triangles and points from file (mesh)
SurfaceDir=CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"VertPenny-300-EqEdges.stl")
(Points,Triangles)=CutAndDisplaceJulia.STLReader(SurfaceDir)

#Flatten so normals point up
X=Points[:,2];
Points[:,2]=Points[:,4];
Points[:,4]=X;

#Pennys angle away from Z. 
Beta=15; 
BetaFromVert=90-Beta;
#Beta=deg2rad(90-Beta)
#Rotate this (XZ)
(Points[:,2],Points[:,4])=CutAndDisplaceJulia.RotateObject2D!(Points[:,2],Points[:,4],0.0,0.0,cosd(BetaFromVert),sind(BetaFromVert))

#Compute triangle properties
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
n=length(Triangles[:,1]);
n2=length(Points[:,1]);

(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )

# Which bits we want to compute
HSFlag=0; #const 

#Elastic constants
G=ShearModulus(5.); 
ν=PoissonsRatio(0.25);
(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(G,ν);

#Set BoundaryConditions
σxx = zeros(n);	
σzz = ones(n);	
σyy = zeros(n);	
σxy = zeros(n);    
σxz = zeros(n);
σyz = zeros(n);
BoundaryConditions=Stresses(σxx,σyy,σzz,σxy,σxz,σyz);

FixedEls=zeros(n,1);

#Calculate slip on faces
(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);

#Get tip elements
(FeP1P2S,FeP1P3S,FeP2P3S)=CutAndDisplaceJulia.GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector)
(FeP1P2S,FeP1P3S,FeP2P3S)=CutAndDisplaceJulia.StressIntensity3D(Dn,Dss,Dds,G,ν,FaceNormalVector,FeP1P2S,FeP1P3S,FeP2P3S);

#Calculate theta [location around the crack
FreeEdMdX=[FeP1P2S.FeMd[FeP1P2S.FreeFlg,1];FeP1P3S.FeMd[FeP1P3S.FreeFlg,1];FeP2P3S.FeMd[FeP2P3S.FreeFlg,1]];
FreeEdMdY=[FeP1P2S.FeMd[FeP1P2S.FreeFlg,2];FeP1P3S.FeMd[FeP1P3S.FreeFlg,2];FeP2P3S.FeMd[FeP2P3S.FreeFlg,2]];
FreeEdMdZ=[FeP1P2S.FeMd[FeP1P2S.FreeFlg,3];FeP1P3S.FeMd[FeP1P3S.FreeFlg,3];FeP2P3S.FeMd[FeP2P3S.FreeFlg,3]];
# Flatten back to calc theta


(FreeEdMdX,FreeEdMdZ) = CutAndDisplaceJulia.RotateObject2D(FreeEdMdX,FreeEdMdZ,0,0,cosd(-BetaFromVert),sind(-BetaFromVert));
(θ,~)=CutAndDisplaceJulia.cart2pol(FreeEdMdX,FreeEdMdY);
θ=rad2deg.(θ);
#%Accumulate from structure into big vectors:
K1=[FeP1P2S.K1[FeP1P2S.FreeFlg];FeP1P3S.K1[FeP1P3S.FreeFlg];FeP2P3S.K1[FeP2P3S.FreeFlg]];
K2=[FeP1P2S.K2[FeP1P2S.FreeFlg];FeP1P3S.K2[FeP1P3S.FreeFlg];FeP2P3S.K2[FeP2P3S.FreeFlg]];
K3=[FeP1P2S.K3[FeP1P2S.FreeFlg];FeP1P3S.K3[FeP1P3S.FreeFlg];FeP2P3S.K3[FeP2P3S.FreeFlg]];

#Compute across each theta
(K1an,K2an,K3an) = CutAndDisplaceJulia.Tada_StrIntInclinedPennyTension(σzz[1],Beta,θ,1,ν);

#Compute the percent error between analytical and numerical
ResidualPercentK1=CutAndDisplaceJulia.BAsPercentOfA(K1an,K1);
ResidualPercentK2=CutAndDisplaceJulia.BAsPercentOfA(K2an,K2);
ResidualPercentK3=CutAndDisplaceJulia.BAsPercentOfA(K3an,K3);

@info ResidualPercentK1 ResidualPercentK2 ResidualPercentK3
MaxErrorK1=maximum(filter(!isnan,abs.(ResidualPercentK1)))-100; 
MaxErrorK2=maximum(filter(!isnan,abs.(ResidualPercentK2)))-100; 
MaxErrorK3=maximum(filter(!isnan,abs.(ResidualPercentK3)))-100;
@info MaxErrorK1 MaxErrorK2 MaxErrorK3
#Test this has not changed 
lim=11; #Percent error limit
if MaxErrorK1>lim || MaxErrorK2>lim || MaxErrorK3>lim
 	error("Residual displacement too high, some over $lim%")
end

println("Test Passed")


#=
#To Draw
using Plots
gr()
y=[K1an K1];
plot1=scatter(θ,y,title="Tht vs KI, An=y1 BEM=y2 NOTE Y-AXIS!")
y=[K2an K2];
plot2=scatter(θ,y,title="Tht vs KII, An=y1 BEM=y2")
y=[K3an K3];
plot3=scatter(θ,y,title="Tht vs KIII, An=y1 BEM=y2")
plot(plot1,plot2,plot3,layout=(3,1))
=#