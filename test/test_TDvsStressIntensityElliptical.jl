#Test case: stress intensity (KI) around an elliptical crack (semi-axes a=2 along x, b=1 along y)
#under uniform tension σzz, against the analytical solution (Tada, Stress analysis handbook 26.2).
#The mesh is made once with CGAL by MakeSIFTestMeshes.jl (remeshed, rim cleaned, rim snapped
#onto the ellipse) and stored in this folder, so this test does not need BuildCGAL.
using Test
using Statistics

#Load triangles and points from file (mesh)
SurfaceDir=CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"EllipseMesh_a2_b1_1500Faces.off")
(Points,Triangles)=CutAndDisplaceJulia.OFFReader(SurfaceDir)
(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3(Triangles,Points)
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)

#Semiaxes (a must be bigger than b for the Tada solution)
a=2;
b=1;

n=length(Triangles[:,1]);
HSFlag=0; #full space

#Elastic constants
G=ShearModulus(5.);
ν=PoissonsRatio(0.25);
(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(G,ν);

#Uniform tension perpendicular to the crack
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

#Stress intensity on the tip elements, then smoothed along the rim (as in the propagation loop)
(FeP1P2S,FeP1P3S,FeP2P3S)=CutAndDisplaceJulia.GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector)
(FeP1P2S,FeP1P3S,FeP2P3S)=CutAndDisplaceJulia.StressIntensity3D(Dn,Dss,Dds,G,ν,FaceNormalVector,FeP1P2S,FeP1P3S,FeP2P3S);
FreeEdMdX=[FeP1P2S.FeMd[FeP1P2S.FreeFlg,1];FeP1P3S.FeMd[FeP1P3S.FreeFlg,1];FeP2P3S.FeMd[FeP2P3S.FreeFlg,1]];
FreeEdMdY=[FeP1P2S.FeMd[FeP1P2S.FreeFlg,2];FeP1P3S.FeMd[FeP1P3S.FreeFlg,2];FeP2P3S.FeMd[FeP2P3S.FreeFlg,2]];
avgeverynth=3
(FeP1P2S,FeP1P3S,FeP2P3S)=CutAndDisplaceJulia.MovingAverageOfStressIntensity(avgeverynth,P1,P2,P3,FaceNormalVector,MidPoint,FeP1P2S,FeP1P3S,FeP2P3S)
K1=[FeP1P2S.K1[FeP1P2S.FreeFlg];FeP1P3S.K1[FeP1P3S.FreeFlg];FeP2P3S.K1[FeP2P3S.FreeFlg]];

#Analytical KI at the same places round the rim
(K1an,θ)=CutAndDisplaceJulia.Tada_StrIntEllipseCrackTension(σzz[1],a,b,FreeEdMdX,FreeEdMdY)

#Compare
RelErr=abs.(K1.-K1an)./maximum(K1an)
println("Elliptical crack in tension, a/b = $a/$b, $n triangles: KI error (fraction of max KI) median $(round(median(RelErr),sigdigits=3)), max $(round(maximum(RelErr),sigdigits=3))")
@test median(RelErr)<0.03     #measured: ~0.01-0.015
@test maximum(RelErr)<0.08     #measured: ~0.05
