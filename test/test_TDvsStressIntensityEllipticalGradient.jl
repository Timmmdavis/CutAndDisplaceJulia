#Test case: stress intensity (KI) around an elliptical crack (semi-axes xc=1 along x, zc=3 along y)
#under internal pressure p0 plus a linear gradient along y, sized by CompP0P1 so that KI = Kc at
#the top tip and KI = 0 at the base tip, against the analytical solution of Atroshchenko et al.
#(2009). For this tall crack KI is NOT largest at the top: the pressure part is largest at the
#sides and the gradient part at the top, so the total peaks at ~1.28 Kc about 45° from the top
#(as in Tim's MATLAB EllipticalCrackAssymGrad_K / CrackOpeningEllipse.m, which also warns that
#rim K exceeds Kc there).
#(src/Atroshchenko_StrIntEllipseCrackGradient.jl).
#The mesh is made once with CGAL by MakeSIFTestMeshes.jl and stored in this folder.
using Test
using Statistics

SurfaceDir=CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"EllipseMesh_a1_b3_650Faces.off")
(Points,Triangles)=CutAndDisplaceJulia.OFFReader(SurfaceDir)
(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3(Triangles,Points)
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)

#Semiaxes
xc=1;
zc=3;

n=length(Triangles[:,1]);
HSFlag=0; #full space

#Elastic constants
G=ShearModulus(5.);
ν=PoissonsRatio(0.25);
(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(G,ν);

#Pressure p0 and gradient p1·zc (P1c) giving KI = Kc at the top and 0 at the base
Kc=1;
(c1,c2,k,k_prime,Kk,Ek,Q)=CutAndDisplaceJulia.Assign_c1c2_andEllipticIntegrals(xc,zc)
(p0,p1)=CutAndDisplaceJulia.CompP0P1(c1,c2,zc,xc,Ek,Kk,k,k_prime,Q,Kc)
p1x=p1/zc
σxx = zeros(n);
σzz = ones(n)*p0.+MidPoint[:,2]*(p1x);
σyy = zeros(n);
σxy = zeros(n);
σxz = zeros(n);
σyz = zeros(n);
BoundaryConditions=Stresses(σxx,σyy,σzz,σxy,σxz,σyz);
FixedEls=zeros(n,1);

#Calculate slip on faces
(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);

#Stress intensity on the tip elements, then smoothed along the rim
(FeP1P2S,FeP1P3S,FeP2P3S)=CutAndDisplaceJulia.GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector)
(FeP1P2S,FeP1P3S,FeP2P3S)=CutAndDisplaceJulia.StressIntensity3D(Dn,Dss,Dds,G,ν,FaceNormalVector,FeP1P2S,FeP1P3S,FeP2P3S);
FreeEdMdX=[FeP1P2S.FeMd[FeP1P2S.FreeFlg,1];FeP1P3S.FeMd[FeP1P3S.FreeFlg,1];FeP2P3S.FeMd[FeP2P3S.FreeFlg,1]];
FreeEdMdY=[FeP1P2S.FeMd[FeP1P2S.FreeFlg,2];FeP1P3S.FeMd[FeP1P3S.FreeFlg,2];FeP2P3S.FeMd[FeP2P3S.FreeFlg,2]];
avgeverynth=3
(FeP1P2S,FeP1P3S,FeP2P3S)=CutAndDisplaceJulia.MovingAverageOfStressIntensity(avgeverynth,P1,P2,P3,FaceNormalVector,MidPoint,FeP1P2S,FeP1P3S,FeP2P3S)
K1=[FeP1P2S.K1[FeP1P2S.FreeFlg];FeP1P3S.K1[FeP1P3S.FreeFlg];FeP2P3S.K1[FeP2P3S.FreeFlg]];

#Analytical KI round the rim. Phi = angle from the top tip, measured on the ellipse scaled to a
#circle (divide y by zc and x by xc first)
(θ,ρ)=CutAndDisplaceJulia.cart2pol(FreeEdMdY./zc,FreeEdMdX./xc)
K1an=zeros(size(θ))
for i=1:length(θ)
	(KI_InternalPressure,KI_GradientPressure)=CutAndDisplaceJulia.Comp_KC_EllipseEdge(c1,c2,zc,xc,Ek,Kk,k,k_prime,Q,p0,p1,θ[i])
	K1an[i]=KI_InternalPressure+KI_GradientPressure
end

#The analytical conditions themselves: Kc at the top, 0 at the base
KAt(ph)=sum(CutAndDisplaceJulia.Comp_KC_EllipseEdge(c1,c2,zc,xc,Ek,Kk,k,k_prime,Q,p0,p1,ph))
@test KAt(0.0) ≈ Kc
@test abs(KAt(π))<1e-12

#Numerical KI at the tip elements nearest the top, side and base
Nearest(ph)=argmin(abs.(mod.(θ.-ph.+π,2π).-π))
(Itop,Iside,Ibase)=(Nearest(0.0),Nearest(π/2),Nearest(π))
println("Elliptical crack with pressure gradient, xc/zc = $xc/$zc, $n triangles (KI/Kc, numerical vs analytical):")
println("  top  $(round(K1[Itop],digits=3)) vs $(round(K1an[Itop],digits=3)),  side $(round(K1[Iside],digits=3)) vs $(round(K1an[Iside],digits=3)),  base $(round(K1[Ibase],digits=3)) vs $(round(K1an[Ibase],digits=3))")
println("  largest: $(round(maximum(K1),digits=3)) at $(round(rad2deg(θ[argmax(K1)]),digits=0))° vs $(round(maximum(K1an),digits=3)) at $(round(rad2deg(θ[argmax(K1an)]),digits=0))° from the top")
@test abs(K1[Itop]-K1an[Itop])<0.12*Kc     #measured ~0.09 (the top is the most curved part of the rim)
@test abs(K1[Ibase])<0.1*Kc                #measured ~0.03
@test abs(abs(rad2deg(θ[argmax(K1)]))-abs(rad2deg(θ[argmax(K1an)])))<15   #peak in the same place

#Whole rim
RelErr=abs.(K1.-K1an)./Kc
println("  error round the rim (fraction of Kc): median $(round(median(RelErr),sigdigits=3)), max $(round(maximum(RelErr),sigdigits=3))")
@test median(RelErr)<0.03     #measured: ~0.015
@test maximum(RelErr)<0.15     #measured: ~0.10 (largest near the peak)
