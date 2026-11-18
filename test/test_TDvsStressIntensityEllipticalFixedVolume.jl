# Fixed-volume (incompressible fluid) elliptical cracks in a stress gradient, against the analytical solution.
#
# The fluid volume solve (SlipCalculator3D with MixedBoundaryConditionsFluidVolume) is given ONLY the stress
# gradient across the crack (the fluid weight, say) and the crack VOLUME; it must find the fluid pressure itself.
# The volume comes from Atroshchenko et al. (2009) (src/Atroshchenko_StrIntEllipseCrackGradient.jl): the uniform
# pressure p0 and gradient that give KI = 0 at the base tip (-y) and KI = Kc at the top tip (+y). The gradient part
# is odd in y and adds no volume, so  V = 8 pi p0 c1 c2^2 (1-nu^2) / (3 E E(k))  (c1, c2 the long and short semi-axes).
#
# Three shapes (semi-axis a along x, b along the gradient y): oblate 2/1, circular 1/1, prolate 1/3. Checks:
#  1. the volume that comes out is the volume that went in (solver precision);
#  2. the pressure the solver found is the analytical p0;
#  3. KI round the rim is ~0 at the base, ~Kc at the top, and follows the analytical curve.
# 2 and 3 hold to the accuracy of the BEM mesh (the tolerances below).
# The meshes are built here with GradedRemesh (CGAL), as the propagation loop's graded meshes: small triangles along
# the rim, growing inwards. The rim is put exactly on the ellipse and kept there (checked). The test needs BuildCGAL
# with its graded_remeshing tool and is skipped without it (as unit test 11).
# Figure: each crack coloured by opening (top row), and KI/Kc round its rim against the analytical curve (bottom).
include(joinpath(@__DIR__,"UnitTests","UnitTestHelpers.jl"))
using Printf
const CDJ=CutAndDisplaceJulia

#Graded ellipse mesh, semi-axes a (x) and b (y), in the x-y plane (normals +z), ~NoTris triangles. The rim is made
#first, exactly on the ellipse: points equally spaced by arc length at the rim edge size hf, joined to the centre by
#a fan. GradedRemesh (CGAL) then meshes the inside, keeping that rim as given: edges hf on the rim growing to h inside,
#and one rim triangle per rim edge (no clean-up needed). Sizes as in the propagation loop's graded scheme: h = 1.4x
#and hf = 0.7x the edge of NoTris equilateral triangles covering the ellipse.
function EllipseMesh(a,b,NoTris)
	h0=sqrt(4π*a*b/(sqrt(3)*NoTris)); h=1.4*h0; hf=0.7*h0
	t=range(0,2π,length=20001)
	L=[0;cumsum(hypot.(diff(a.*cos.(t)),diff(b.*sin.(t))))]       #arc length along the rim
	N=round(Int,L[end]/hf)
	tk=[t[searchsortedfirst(L,l)] for l in range(0,L[end],length=N+1)[1:N]]
	Points=[collect(1:N+1) [0.0;a.*cos.(tk)] [0.0;b.*sin.(tk)] zeros(N+1)]
	Triangles=[fill(1,N) collect(2:N+1) [collect(3:N+1);2]]
	(Points,Triangles)=cd(mktempdir()) do                        #the CGAL round trip writes .off files
		CDJ.GradedRemesh(Points,Triangles,h,hf,BuildCGAL.GradedRemeshingCGAL;growth=1.3)
	end
	(P1,P2,P3)=CDJ.CreateP1P2P3(Triangles,Points)
	(FaceNormalVector,MidPoint)=CDJ.CreateFaceNormalAndMidPoint(Points,Triangles)
	return P1,P2,P3,FaceNormalVector,MidPoint
end

#Largest distance of the mesh's rim points off the ellipse, as |x²/a² + y²/b² - 1|
function RimOffEllipse(P1,P2,P3,a,b)
	E=BoundaryEdgesBruteForce(P1,P2,P3)
	return maximum(abs((p[1]/a)^2+(p[2]/b)^2-1) for e in E for p in e)
end

#Analytical solution: p0 and the gradient (stress per unit y) for KI = 0 at the base and Kc at the top,
#the crack volume, and KI at angle θ round the rim (from the top tip)
function EllipseAnalytical(a,b,Kc,E,ν)
	(c1,c2,k,k_prime,Kk,Ek,Q)=CDJ.Assign_c1c2_andEllipticIntegrals(a,b)
	(p0,P1c)=CDJ.CompP0P1(c1,c2,b,a,Ek,Kk,k,k_prime,Q,Kc)
	KI(θ)=sum(CDJ.Comp_KC_EllipseEdge(c1,c2,b,a,Ek,Kk,k,k_prime,Q,p0,P1c,θ))
	return (p0=p0,Gradient=P1c/b,Volume=8*pi*p0*c1*c2^2*(1-ν^2)/(3*E*Ek),KI=KI)
end

#Numerical KI round the rim (smoothed as in the other KI tests), sorted by angle θ from the top tip, measured on
#the ellipse (0 top, +90° at +x, ±180° base)
function RimKI(Dn,Dss,Dds,P1,P2,P3,FaceNormalVector,MidPoint,G,ν,a,b)
	(E12,E13,E23)=CDJ.GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector)
	(E12,E13,E23)=CDJ.StressIntensity3D(Dn,Dss,Dds,G,ν,FaceNormalVector,E12,E13,E23)
	Rim(c)=[E12.FeMd[E12.FreeFlg,c];E13.FeMd[E13.FreeFlg,c];E23.FeMd[E23.FreeFlg,c]]
	(X,Y)=(Rim(1),Rim(2))
	(E12,E13,E23)=CDJ.MovingAverageOfStressIntensity(3,P1,P2,P3,FaceNormalVector,MidPoint,E12,E13,E23)
	K1=[E12.K1[E12.FreeFlg];E13.K1[E13.FreeFlg];E23.K1[E23.FreeFlg]]
	(θ,~)=CDJ.cart2pol(Y./b,X./a)
	o=sortperm(θ)
	return θ[o],K1[o]
end

#One crack: mesh, the fixed-volume solve with the analytical volume, and its result
function FixedVolumeCase(a,b,NoTris,G,ν,λ,E,Kc)
	(P1,P2,P3,FaceNormalVector,MidPoint)=EllipseMesh(a,b,NoTris)
	n=size(P1,1); Zero=zeros(n); HSFlag=0
	(Area,~)=TriAreas(P1,P2,P3)
	An=EllipseAnalytical(a,b,Kc,E,ν)
	σGrad=An.Gradient.*MidPoint[:,2]
	Plain(σzz)=CDJ.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,Stresses(Zero,Zero,σzz,Zero,Zero,Zero),zeros(n,1))
	FixedVolume(V)=CDJ.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,
		MixedBoundaryConditionsFluidVolume(MixedBoundaryConditions(Stresses(Zero,Zero,σGrad,Zero,Zero,Zero),Tractions(Zero,Zero,Zero)),V),ones(n))
	Rim(D)=RimKI(D...,P1,P2,P3,FaceNormalVector,MidPoint,G,ν,a,b)

	D=FixedVolume(An.Volume)
	(θ,K1)=Rim(D)
	#The pressure the solve found: Dn is linear in it, Dn = P*Dn(unit pressure) + Dn(gradient alone)
	DnUnit=Plain(ones(n))[1]; DnGrad=Plain(σGrad)[1]
	P=sum(DnUnit.*(D[1].-DnGrad))/sum(DnUnit.^2)
	return (a=a,b=b,n=n,P1=P1,P2=P2,P3=P3,Dn=D[1],An=An,θ=θ,K1=K1,K1an=An.KI.(θ),
		Vin=An.Volume,Vout=sum(D[1].*Area),P=P,RimOff=RimOffEllipse(P1,P2,P3,a,b))
end

HasGradedTool=isdefined(Main,:BuildCGAL) && isdefined(BuildCGAL,:GradedRemeshingCGAL) &&
	isfile(joinpath(BuildCGAL.tools_dir(),Sys.iswindows() ? "graded_remeshing.exe" : "graded_remeshing"))
if !HasGradedTool
	@info "BuildCGAL with the graded_remeshing tool not loaded: skipping the fixed-volume elliptical KI test (it builds its meshes with CGAL)"
else
@testset "Elliptical cracks, fixed volume in a stress gradient: KI = 0 at the base, Kc at the top" begin

	G=ShearModulus(5.0); ν=PoissonsRatio(0.25)
	(K,E,λ,ν,G)=CDJ.ElasticConstantsCheck(G,ν)
	Kc=1.0
	NoTris=600
	Cases=[("Oblate",2.0,1.0),("Circular",1.0,1.0),("Prolate",1.0,3.0)]   #(name, a along x, b along y)

	#Tolerances (fractions of p0 / of Kc): what these ~670-triangle graded meshes give, with margin
	TolP=0.07            #|P/p0-1|: measured 0.031, 0.020, 0.043
	TolBase=0.05         #|KI(base)|/Kc: measured 0.000, 0.000, 0.000
	TolTop=0.08          #|KI(top)-Kc|/Kc: measured 0.038, 0.031, 0.025
	TolRim=0.06          #median |KI - KI analytical|/Kc round the rim: measured 0.005, 0.014, 0.044

	Results=[]
	for (Name,a,b) in Cases
		R=FixedVolumeCase(a,b,NoTris,G,ν,λ,E,Kc)
		Itop=argmin(abs.(R.θ)); Ibase=argmax(abs.(R.θ))
		RimErr=median(abs.(R.K1.-R.K1an))/Kc
		push!(Results,merge(R,(Name=Name,Itop=Itop,Ibase=Ibase,RimErr=RimErr)))
		@printf("  %-9s a/b=%g/%g %4d tris, rim off ellipse %.0e | V out/in %.10f | P/p0 %.3f | KI/Kc top %.3f base %+.3f, rim median error %.3f\n",
			Name,a,b,R.n,R.RimOff,R.Vout/R.Vin,R.P/R.An.p0,R.K1[Itop],R.K1[Ibase],RimErr)

		@test R.An.KI(0.0)≈Kc                                  #the analytical loading is what we think it is
		@test abs(R.An.KI(π))<1e-12*Kc
		@test R.RimOff<1e-8                                     #the mesh rim is the ellipse
		@test R.Vout≈R.Vin rtol=1e-8                            #1. volume held
		@test abs(R.P/R.An.p0-1)<TolP                           #2. pressure found
		@test abs(R.K1[Ibase])<TolBase*Kc                       #3. KI against the analytical curve
		@test abs(R.K1[Itop]-Kc)<TolTop*Kc
		@test RimErr<TolRim
		GC.gc()
	end

	if MakeFigures
		fig=CM.Figure(size=(1500,1000))
		L=3.3                                               #same scale for all three cracks
		local plt
		for (i,R) in enumerate(Results)
			ax=CM.Axis(fig[1,i],aspect=CM.DataAspect(),xlabel="x",ylabel=i==1 ? "y (gradient direction)" : "",
				title="$(R.Name): a/b = $(R.a)/$(R.b), $(R.n) triangles")
			plt=FillTris!(ax,R.P1,R.P2,R.P3,R.Dn./maximum(R.Dn);colormap=:viridis,colorrange=(0,1))
			DrawTriEdges!(ax,R.P1,R.P2,R.P3;color=(:white,0.35),linewidth=0.4)
			t=range(0,2π,length=200)
			CM.lines!(ax,R.a.*sin.(t),R.b.*cos.(t),color=:black,linewidth=1)
			CM.scatter!(ax,[0.0],[R.b],color=:red,markersize=14,strokecolor=:black,strokewidth=1)
			CM.scatter!(ax,[0.0],[-R.b],color=:dodgerblue,markersize=14,strokecolor=:black,strokewidth=1)
			CM.text!(ax,0.0,R.b,text="  top: KI = Kc",align=(:left,:bottom),fontsize=13)
			CM.text!(ax,0.0,-R.b,text="  base: KI = 0",align=(:left,:top),fontsize=13)
			CM.limits!(ax,-L,L,-L,L)

			ax2=CM.Axis(fig[2,i],xlabel="angle round the rim from the top tip (°)",ylabel=i==1 ? "KI / Kc" : "",
				xticks=-180:90:180,
				title=@sprintf("V out / V in = %.8f,  pressure found / p0 = %.3f\nKI/Kc: top %.2f (1), base %+.2f (0), rim median error %.3f",
					R.Vout/R.Vin,R.P/R.An.p0,R.K1[R.Itop],R.K1[R.Ibase],R.RimErr),titlesize=13)
			d=rad2deg.(R.θ)
			CM.band!(ax2,d,R.K1an.-0.05,R.K1an.+0.05,color=(:gray,0.25),label="analytical ± 0.05 Kc")
			CM.lines!(ax2,d,R.K1an,color=:black,linewidth=2,label="analytical (Atroshchenko et al. 2009)")
			CM.scatter!(ax2,d,R.K1,color=:blue,markersize=6,label="fixed-volume solve (BEM)")
			CM.scatter!(ax2,[0.0,180.0],[1.0,0.0],color=[:red,:dodgerblue],markersize=14,strokecolor=:black,strokewidth=1)
			CM.limits!(ax2,-185,185,-0.25,1.35)
			i==1 && CM.axislegend(ax2,position=:cb,framevisible=false,labelsize=12)
		end
		CM.Colorbar(fig[1,4],plt,label="opening / maximum opening")
		CM.Label(fig[0,:],"Fixed-volume elliptical cracks in a stress gradient: given only the gradient and the analytical volume,\n"*
			"the solve must find the pressure that gives KI = 0 at the base and Kc at the top",fontsize=19)
		SaveTestFigure(fig,"test_TDvsStressIntensityEllipticalFixedVolume")
	end
end
end
