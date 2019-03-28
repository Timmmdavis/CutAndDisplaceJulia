#Test case comparing to Penny shaped crack

#Start creating vars for function: 
println("creating func vars")

#Load triangles and points from file (mesh)
SurfaceDir=CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"CircleMesh_1a_500Faces.ts")
(Points,Triangles)=CutAndDisplaceJulia.GoCadAsciiReader(SurfaceDir)

#Compute triangle properties
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
n=length(Triangles[:,1]);
n2=length(Points[:,1]);

(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )

#Elastic constants
G=ShearModulus(1.); 
ν=PoissonsRatio(0.25);
(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(G,ν);


# Which bits we want to compute
HSFlag=1; #const 


#Traction vector
Tn=zeros(n,1);
Tds=zeros(n,1);
Tss=ones(n,1);

FixedEls=zeros(n,1);

Tn=0.0;
Tds=0.0;
Tss=1.0;

#Set BoundaryConditions
BoundaryConditions=Tractions(Tn,Tss,Tds)

#Calculate slip on faces
(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);

#Compute total shearing
TotalShearing = sqrt.((Dss).^2 .+(Dds).^2);

#Finding distance from 0,0 for code midpoints for plotting. 
(θ,ρ) = CutAndDisplaceJulia.cart2pol(MidPoint[:,1],MidPoint[:,2]);

#Compute analytical solution for shearing of cracks walls
a=1.0; #Radius of penny
(us)=CutAndDisplaceJulia.Eshelby1957_PennyCrackSlipProfile(G,ν,0.0,1.0,ρ,a)#Shearing
(ut)=CutAndDisplaceJulia.Eshelby1957_PennyCrackSlipProfile(G,ν,1.0,0.0,ρ,a)#Opening

#Compute the percent error between analytical and numerical
ResidualPercent=CutAndDisplaceJulia.BAsPercentOfA(us,TotalShearing);

#Test this has not changed 
using Statistics
if maximum(ResidualPercent)>0.012 || Statistics.mean(ResidualPercent)>0.0055
	error("Residual is too high")
end

println("Test Passed")