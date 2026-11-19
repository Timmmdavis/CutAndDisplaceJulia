# Unit test 13: rock weight + fluid weight give the right buoyancy loading, and the crack opens at the TOP.
# Vertical penny (stored mesh, no remeshing) of the critical radius, centre 5a below the surface (z up).
#  - Loading: rock stress (σzz = ρrock*g*z, σxx = σyy = StressRatio*σzz, tension positive) plus the weight of
#    the fluid as a traction (Tn = -ρfluid*g*z, opening-positive). Resolved on the crack by
#    SetupTractionVector this must equal (StressRatio*ρrock - ρfluid)*g*z, i.e. Δρ*z when isotropic.
#  - Direction: with a normal load that increases upwards (Δρ*(z-zc) + p0) the opening is centred above the
#    middle of the crack and the stress intensity is larger at the top tip than the bottom tip.
#  - Fluid volume: the volume solve (as used by testProp) with the critical volume gives that volume, opens
#    towards the top, and the top tip is near critical (K1 ~ Kc), which checks the critical radius
#    a = (3*sqrt(π)*Kc/(8*Δρ))^(2/3) and volume formulae used by testProp.
include(joinpath(@__DIR__,"UnitTestHelpers.jl"))

@testset "Buoyancy loading: rock weight + fluid weight" begin

	ρrock=2900.0; ρfluid=2600.0; g=9.81
	Δρ=(ρrock-ρfluid)*g                          #[Pa/m]
	Kc=5e7                                       #[Pa√m]
	G=ShearModulus(2.0e9); ν=PoissonsRatio(0.25)
	(K,E,λ,ν,G)=CutAndDisplaceJulia.ElasticConstantsCheck(G,ν)
	a=((3*sqrt(pi)*Kc)/(8*Δρ))^(2/3)             #critical radius
	Vcrit=real(((1-ν)/(16*G))*((9*π^4*Kc^8)/(complex(Δρ)^5))^(1/3))   #critical volume (vertical, isotropic)

	#Vertical penny (normals horizontal), radius a, centre 5a below the surface
	SurfaceDir=CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"VertPenny-300-EqEdges.stl")
	(Points,Triangles)=CutAndDisplaceJulia.STLReader(SurfaceDir)
	Points[:,2:4]=Points[:,2:4].*a
	Points[:,4]=Points[:,4].-5*a
	(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3(Triangles,Points)
	(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
	n=size(P1,1); Z=MidPoint[:,3]; Zc=mean(Z)
	(Area,HalfPerimeter)=TriAreas(P1,P2,P3)
	@test maximum(abs.(FaceNormalVector[:,3]))<1e-6      #normals horizontal
	@test all(Z.<0)                                      #below the surface

	#Rock stress + fluid weight, as in testProp
	#Normalise=true shifts the net normal load so its maximum is zero, as testProp does before the volume solve
	#(the load is about -Δρ*depth, far beyond the pressure range the volume solve searches; a uniform shift is
	#absorbed by the uniform pressure it solves for)
	function BuoyancyLoading(StressRatio; Normalise=false)
		σzz=ρrock*g.*Z
		σh=StressRatio.*σzz
		Stress=Stresses(σh,σh,σzz,zeros(n),zeros(n),zeros(n))
		TnFluid=-ρfluid*g.*Z
		if Normalise
			(TnStress,~,~)=CutAndDisplaceJulia.CalculateNormalAndShearTractions3D(σh,σh,σzz,zeros(n),zeros(n),zeros(n),FaceNormalVector)
			TnFluid=TnFluid.-maximum(TnFluid.+TnStress)
		end
		return MixedBoundaryConditions(Stress,Tractions(TnFluid,zeros(n),zeros(n)))
	end

	#Stress intensity at the tips: tip positions and K1
	function TipStressIntensity(Dn,Dss,Dds)
		(F1,F2,F3)=CutAndDisplaceJulia.GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector)
		(F1,F2,F3)=CutAndDisplaceJulia.StressIntensity3D(Dn,Dss,Dds,G,ν,FaceNormalVector,F1,F2,F3)
		Fe=(F1,F2,F3)
		return vcat([F.FeMd[F.FreeFlg,:] for F in Fe]...), vcat([F.K1[F.FreeFlg] for F in Fe]...)
	end
	Centroid(Dn)=sum(Area.*Dn.*Z)/sum(Area.*Dn)          #height of the middle of the opening

	#1. The net load on a vertical crack
	for StressRatio in (1.0,0.6)
		(TnNet,~,~)=CutAndDisplaceJulia.SetupTractionVector(BuoyancyLoading(StressRatio),FaceNormalVector,n,G,λ)
		@test TnNet ≈ (StressRatio*ρrock-ρfluid)*g.*Z rtol=1e-6
	end
	(TnIso,~,~)=CutAndDisplaceJulia.SetupTractionVector(BuoyancyLoading(1.0),FaceNormalVector,n,G,λ)
	@test TnIso ≈ Δρ.*Z rtol=1e-6                         #isotropic: buoyancy, larger (less negative) upwards

	#2. A load that increases upwards opens the crack most at the top
	p0=1.5*Δρ*a                                            #keeps the load positive (open) everywhere
	TnUp=Δρ.*(Z.-Zc).+p0
	(Dn,Dss,Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,0,Tractions(TnUp,zeros(n),zeros(n)),zeros(n))
	@test all(Dn.>0)
	@test Centroid(Dn)>Zc
	(TipXYZ,KI)=TipStressIntensity(Dn,Dss,Dds)
	Top=TipXYZ[:,3].>Zc+0.5*a; Bottom=TipXYZ[:,3].<Zc-0.5*a
	@test mean(KI[Top])>mean(KI[Bottom])
	#...and the opposite load (increasing downwards) opens it most at the bottom: the test can tell them apart
	(DnDown,DssDown,DdsDown)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,0,Tractions(-Δρ.*(Z.-Zc).+p0,zeros(n),zeros(n)),zeros(n))
	@test Centroid(DnDown)<Zc

	#3. Fixed fluid volume (as testProp): the critical volume, rock + fluid weight
	BC=MixedBoundaryConditionsFluidVolume(BuoyancyLoading(1.0;Normalise=true),Vcrit)
	(Dn,Dss,Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,0,BC,fill(1,n))
	@test sum(Area.*Dn) ≈ Vcrit rtol=0.02
	@test Centroid(Dn)>Zc
	(TipXYZ,KI)=TipStressIntensity(Dn,Dss,Dds)
	Top=TipXYZ[:,3].>Zc+0.5*a; Bottom=TipXYZ[:,3].<Zc-0.5*a
	KTop=maximum(KI[Top]); KBottom=maximum(KI[Bottom])
	@info "Critical radius a, volume Vcrit; K1 top and bottom tip over Kc" a Vcrit KTop/Kc KBottom/Kc
	@test KTop>KBottom
	@test 0.5<KTop/Kc<2                                    #the top tip is close to critical (checks a and Vcrit)
end
