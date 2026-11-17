# Unit test 5: stress tensors and tractions on planes.
#  - CalculateNormalAndShearTractions3D: tractions on planes of every dip under σxx = -1, σzz = -3
#    (tension positive) match the textbook 2D result (Mohr circle):
#      Tn = σxx sin²δ + σzz cos²δ,   |Ts| = |σxx - σzz| sinδ cosδ      (δ = dip of the plane)
#  - TensorTransformation3D: a rotated tensor keeps its invariants and equals A·σ·A'
#  - EigCalc3D: principal values equal the eigenvalues, and σ·dir = S·dir
#  - HookesLaw3DStrain2Stress: uniaxial strain and pure shear give the textbook stresses
#  - cart2pol / pol2cart: round trip
# Figure: Mohr circle with the computed tractions on it, and Tn, |Ts| against dip.
include(joinpath(@__DIR__,"UnitTestHelpers.jl"))

@testset "Stress and tractions" begin

	#Tractions on planes dipping δ = 0..180°, normals in the x-z plane
	σxx=-1.0; σzz=-3.0
	δ=collect(0.0:5.0:180.0)
	n=length(δ)
	Normals=[sind.(δ) zeros(n) cosd.(δ)]
	(Tn,Tds,Tss)=CutAndDisplaceJulia.CalculateNormalAndShearTractions3D(fill(σxx,n),zeros(n),fill(σzz,n),zeros(n),zeros(n),zeros(n),Normals)
	Ts=sqrt.(Tds.^2 .+Tss.^2)
	TnAn=σxx.*sind.(δ).^2 .+σzz.*cosd.(δ).^2
	TsAn=abs.((σxx-σzz).*sind.(δ).*cosd.(δ))
	@test maximum(abs.(Tn.-TnAn))<1e-12
	@test maximum(abs.(Ts.-TsAn))<1e-12

	#A general symmetric tensor and a rotation (30° about z then 20° about x)
	σ=[2.0 0.5 -0.3; 0.5 -1.0 0.8; -0.3 0.8 0.4]
	Rz=[cosd(30) -sind(30) 0; sind(30) cosd(30) 0; 0 0 1]
	Rx=[1 0 0; 0 cosd(20) -sind(20); 0 sind(20) cosd(20)]
	A=Rx*Rz
	(Txx,Tyy,Tzz,Txy,Txz,Tyz)=CutAndDisplaceJulia.TensorTransformation3D(σ[1,1],σ[2,2],σ[3,3],σ[1,2],σ[1,3],σ[2,3],A)
	σ2=[Txx Txy Txz; Txy Tyy Tyz; Txz Tyz Tzz]
	@test tr(σ2) ≈ tr(σ)
	@test det(σ2) ≈ det(σ)
	@test σ2 ≈ A*σ*A' atol=1e-12

	#Principal stresses
	(S1,S2,S3,S1dir,S2dir,S3dir)=CutAndDisplaceJulia.EigCalc3D([σ[1,1]],[σ[2,2]],[σ[3,3]],[σ[1,2]],[σ[1,3]],[σ[2,3]])
	S=[S1[1],S2[1],S3[1]]
	@test sort(S) ≈ sort(eigvals(Symmetric(σ)))
	for (s,d) in zip(S,(S1dir,S2dir,S3dir))
		d=vec(d)
		@test σ*d ≈ s*d atol=1e-10
	end

	#Hooke's law
	G=1.0; ν=0.25
	(K,E,λ,ν,G)=CutAndDisplaceJulia.ElasticConstantsCheck(ShearModulus(G),PoissonsRatio(ν))
	ε=1e-3
	(Sxx,Syy,Szz,Sxy,Sxz,Syz)=CutAndDisplaceJulia.HookesLaw3DStrain2Stress(ε,0.0,0.0,0.0,0.0,0.0,λ,G)
	@test Sxx ≈ (λ+2G)*ε
	@test Syy ≈ λ*ε && Szz ≈ λ*ε
	(Sxx,Syy,Szz,Sxy,Sxz,Syz)=CutAndDisplaceJulia.HookesLaw3DStrain2Stress(0.0,0.0,0.0,ε,0.0,0.0,λ,G)
	@test Sxy ≈ 2G*ε
	@test abs(Sxx)+abs(Syy)+abs(Szz)<1e-15

	#Polar coordinates
	x=[1.0,-2.0,0.5]; y=[0.0,1.0,-3.0]
	(θ,ρ)=CutAndDisplaceJulia.cart2pol(x,y)
	(x2,y2)=CutAndDisplaceJulia.pol2cart(θ,ρ)
	@test x2 ≈ x && y2 ≈ y
	@test ρ ≈ hypot.(x,y)

	if MakeFigures
		fig=CM.Figure(size=(1300,560))
		ax1=CM.Axis(fig[1,1],aspect=CM.DataAspect(),xlabel="normal traction Tn (tension +)",ylabel="shear traction |Ts|",
			title="Mohr circle for σxx = $σxx, σzz = $σzz\nline: analytical, dots: CalculateNormalAndShearTractions3D\n(dips δ and 180°−δ land on the same point)")
		c=(σxx+σzz)/2; r=abs(σxx-σzz)/2
		t=range(0,π,length=200)
		CM.lines!(ax1,c.+r.*cos.(t),r.*sin.(t),color=:black)
		sc=CM.scatter!(ax1,Tn,Ts,color=δ,colormap=:viridis,markersize=10)
		CM.Colorbar(fig[1,2],sc,label="dip of plane δ [°]")
		ax2=CM.Axis(fig[1,3],xlabel="dip of plane δ [°] (normal = (sin δ, 0, cos δ))",ylabel="traction",
			title="Tractions against dip\nmax error: Tn $(round(maximum(abs.(Tn.-TnAn)),sigdigits=2)), |Ts| $(round(maximum(abs.(Ts.-TsAn)),sigdigits=2))",
			xticks=0:45:180)
		CM.lines!(ax2,δ,TnAn,color=:blue,label="Tn analytical")
		CM.scatter!(ax2,δ,Tn,color=:blue,markersize=8,label="Tn code")
		CM.lines!(ax2,δ,TsAn,color=:red,label="|Ts| analytical")
		CM.scatter!(ax2,δ,Ts,color=:red,markersize=8,label="|Ts| code")
		CM.axislegend(ax2,position=:rc)
		CM.Label(fig[2,:],"Also checked (no figure): TensorTransformation3D keeps trace/determinant and equals AσA'; EigCalc3D principal values/directions; Hooke's law uniaxial strain and pure shear; cart2pol/pol2cart round trip",
			fontsize=12,tellwidth=false)
		CM.Label(fig[0,:],"Unit test 5: stresses, principal stresses and tractions on planes",fontsize=20)
		SaveTestFigure(fig,"ut_05_StressAndTractions")
	end
end
