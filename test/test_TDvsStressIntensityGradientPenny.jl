#Test case: stress intensity (KI) around a vertical penny-shaped crack opened by a buoyant fluid:
#internal pressure P plus a gradient Δρ·z, sized so the top tip is critical (KI = KCrit) and the
#base is at KI ≈ 0. Analytical: Tada uniform pressure + Tada linear gradient (penny).
#The mesh (unit radius, x-z plane) is made once with CGAL by MakeSIFTestMeshes.jl and stored in
#this folder; it is scaled to the critical radius here.
using Test
using Statistics

#Elastic constants
G=ShearModulus(2.0e9);
ν=PoissonsRatio(0.25);
(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(G,ν);

#Density
ρrock=2900;
ρfluid=2600;
g=9.81;
Δρ=((ρrock-ρfluid)*g)
KCrit=5e7; #[5e7 = 50 MPa √m]

#Radius at which the top tip is critically stressed, and the pressure at the centre
CritRadius=abs(real((complex(9*Δρ)^(1/3)*(pi*(KCrit^2))^(1/3)))/(4*Δρ))
P=(2*Δρ*CritRadius)/3

#Load the unit vertical penny and scale it to CritRadius
SurfaceDir=CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"VertPennyMesh_r1_1500Faces.off")
(Points,Triangles)=CutAndDisplaceJulia.OFFReader(SurfaceDir)
Points[:,2:4]=Points[:,2:4].*CritRadius
a=CritRadius;
b=CritRadius;
(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3(Triangles,Points)
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)

n=length(Triangles[:,1]);
HSFlag=0; #full space
Z=MidPoint[:,3];

#Boundary conditions: no remote stress, fluid pressure P + Δρ·z on the walls (buoyancy)
σxx = zeros(n);
σyy = zeros(n);
σzz = zeros(n);
σxy = zeros(n);
σxz = zeros(n);
σyz = zeros(n);
Tn=Δρ.*Z.+P
Tss=zeros(n)
Tds=zeros(n)
Stress=Stresses(σxx,σyy,σzz,σxy,σxz,σyz);
Traction=Tractions(Tn,Tss,Tds);
BoundaryConditions=MixedBoundaryConditions(Stress,Traction)
FixedEls=zeros(n,1);

#Calculate slip on faces
(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);

#Stress intensity on the tip elements, then smoothed along the rim
(FeP1P2S,FeP1P3S,FeP2P3S)=CutAndDisplaceJulia.GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector)
(FeP1P2S,FeP1P3S,FeP2P3S)=CutAndDisplaceJulia.StressIntensity3D(Dn,Dss,Dds,G,ν,FaceNormalVector,FeP1P2S,FeP1P3S,FeP2P3S);
FreeEdMdX=[FeP1P2S.FeMd[FeP1P2S.FreeFlg,1];FeP1P3S.FeMd[FeP1P3S.FreeFlg,1];FeP2P3S.FeMd[FeP2P3S.FreeFlg,1]];
FreeEdMdY=[FeP1P2S.FeMd[FeP1P2S.FreeFlg,2];FeP1P3S.FeMd[FeP1P3S.FreeFlg,2];FeP2P3S.FeMd[FeP2P3S.FreeFlg,2]];
FreeEdMdZ=[FeP1P2S.FeMd[FeP1P2S.FreeFlg,3];FeP1P3S.FeMd[FeP1P3S.FreeFlg,3];FeP2P3S.FeMd[FeP2P3S.FreeFlg,3]];
θnum=atan.(FreeEdMdX,FreeEdMdZ);
avgeverynth=3
(FeP1P2S,FeP1P3S,FeP2P3S)=CutAndDisplaceJulia.MovingAverageOfStressIntensity(avgeverynth,P1,P2,P3,FaceNormalVector,MidPoint,FeP1P2S,FeP1P3S,FeP2P3S)
K1=[FeP1P2S.K1[FeP1P2S.FreeFlg];FeP1P3S.K1[FeP1P3S.FreeFlg];FeP2P3S.K1[FeP2P3S.FreeFlg]];

#Analytical: uniform pressure part + gradient part
(K1an,θ)=CutAndDisplaceJulia.Tada_StrIntEllipseCrackTension(P,a,b,FreeEdMdX,FreeEdMdY)
K1anx=CutAndDisplaceJulia.Tada_StrIntPennyGradient(Δρ,a,θnum)
K1an=K1an.+K1anx;

#Compare (KI runs from 0 at the base to KCrit at the top)
RelErr=abs.(K1.-K1an)./KCrit
println("Vertical buoyant penny, r = $(round(CritRadius,digits=1)) m, $n triangles: KI error (fraction of KCrit) median $(round(median(RelErr),sigdigits=3)), max $(round(maximum(RelErr),sigdigits=3))")
println("  top: numerical KI/KCrit $(round(maximum(K1)/KCrit,sigdigits=3)), analytical $(round(maximum(K1an)/KCrit,sigdigits=3))")
@test median(RelErr)<0.03     #measured: ~0.01-0.015
@test maximum(RelErr)<0.08     #measured: ~0.05
