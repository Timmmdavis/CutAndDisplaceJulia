#Test case comparing to Penny shaped crack

#Start creating vars for function: 
println("creating func vars")

#Load triangles and points from file (mesh)
SurfaceDir=CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"CircleMesh_1a_500Faces.ts")
(Points,Triangles)=CutAndDisplaceJulia.GoCadAsciiReader(SurfaceDir)

##TestWithSlightDip
#Points[:,4]=Points[:,2].*1e-9;

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
HSFlag=0; #const 


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


@info ResidualPercent

#Only check values that are not within 10% of the tipline
Good=ρ.<0.9;
ResidualPercent=ResidualPercent[Good];
@info ResidualPercent

#Test this has not changed 
lim=10; #Percent error limit
if any((abs.(ResidualPercent.-100)).>lim) 
 	error("Residual displacement too high, some over $lim%")
end

println("Test Passed")

#To Draw
#using Plots
#gr()
#y=[us TotalShearing];
#scatter(ρ,y,title="R vs Disp, An=y1 BEM=y2")