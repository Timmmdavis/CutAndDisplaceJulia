# Unit test 12: the fluid volume solve (FluidPressureSolver.jl).
#  End to end, through SlipCalculator3D on the small penny:
#   A) uniform compression 1 and the volume V1 that unit tension gives: the crack is fully open and the problem is
#      linear, so the fluid pressure must be 2 and the opening exactly that of unit tension (DnRef).
#   B) compression increasing downwards, 0.3 V1 and B2) 0.6 V1: part of the crack closes; the volume must be held
#      and the opening must not be negative.
#   C) two cracks in one solve (two pressures, two volume rows): each crack must hold its own volume.
#  Solver level, on a dense LCP with the same block structure as the friction LCP (Kaven 2012) and several
#  load cases that close part of the crack: the default solve (reduced system) and the full system must
#  both converge to the same x and pressure, which satisfy the complementarity conditions and the volume.
# Figure: case A opening against unit tension, and case B2 opening.
include(joinpath(@__DIR__,"UnitTestHelpers.jl"))
using Random

#Dense LCP with the block structure of the friction LCP: smooth SPD compliance C (elements on a disk),
#[a] as in SlipCalculator3D(Scl,n,InvertedInfMatA,b,µ,Sf) with µ=0.1, Sf=0, normal load Tn0+grad*x.
function FluidTestLCP(n,seed,Tn0,grad,Vfrac)
	rng=MersenneTwister(seed)
	r=sqrt.(rand(rng,n)); θ=2π*rand(rng,n); X=[r.*cos.(θ) r.*sin.(θ)]
	K=[exp(-sum(abs2,X[i,:]-X[j,:])/0.05) for i=1:n, j=1:n]+0.05I
	C=zeros(3n,3n)
	for k=0:2
		C[k*n+1:(k+1)*n,k*n+1:(k+1)*n]=(k==0 ? 1.0 : 0.8)*K
	end
	C[1:n,n+1:2n]=0.05*K; C[n+1:2n,1:n]=0.05*K
	L1=1:n; L2=n+1:2n; L3=2n+1:3n
	dµ=Diagonal(fill(0.1,n)); ID=Matrix(1.0I,n,n); ZE=zeros(n,n)
	a=[(C[L1,L1]-C[L1,L2]*dµ-C[L1,L3]*dµ) C[L1,L2] C[L1,L3] ZE ZE;
	   (C[L2,L1]-C[L2,L2]*dµ-C[L2,L3]*dµ) C[L2,L2] C[L2,L3] ID ZE;
	   (C[L3,L1]-C[L3,L2]*dµ-C[L3,L3]*dµ) C[L3,L2] C[L3,L3] ZE ID;
	   2*Matrix(dµ) -ID ZE ZE ZE; 2*Matrix(dµ) ZE -ID ZE ZE]
	t=[Tn0 .+ grad*X[:,1]; 0.5*(1 .+ X[:,2]); 0.15*randn(rng,n)]
	b=[C*t; zeros(2n)]
	g=reshape([C*[ones(n);zeros(2n)]; zeros(2n)],:,1)
	w=zeros(5n,1); w[1:n,1].=1/n
	V=[Vfrac*2*(w'*g)[1]]
	return a,b,g,w,V
end

@testset "Fluid volume solver" begin

	(Points,Triangles,P1,P2,P3,FaceNormalVector,MidPoint)=LoadSmallPennyMesh()
	n=size(P1,1)
	(Area,~)=TriAreas(P1,P2,P3)
	G=ShearModulus(1.0); ν=PoissonsRatio(0.25)
	(K,E,λ,ν,G)=CutAndDisplaceJulia.ElasticConstantsCheck(G,ν)
	HSFlag=0
	Zero=zeros(n)
	FluidBC(σzz,Volume)=MixedBoundaryConditionsFluidVolume(
		MixedBoundaryConditions(Stresses(Zero,Zero,σzz,Zero,Zero,Zero),Tractions(Zero,Zero,Zero)),Volume)
	(DnRef,~,~)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,Tractions(ones(n),Zero,Zero),zeros(n))
	V1=sum(DnRef.*Area)
	y=MidPoint[:,2]

	#A: fully open, linear
	(DnA,DssA,DdsA)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,FluidBC(fill(-1.0,n),V1),ones(n))
	@test maximum(abs.(DnA.-DnRef))<1e-6*maximum(DnRef)

	#B, B2: partly closed
	DnB=Dict{String,Vector{Float64}}()
	for (name,Vf) in (("B",0.3),("B2",0.6))
		(Dn,Dss,Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,FluidBC(-1.0.+2.0.*y,Vf*V1),ones(n))
		@test abs(sum(Dn.*Area)/(Vf*V1)-1)<1e-8
		@test minimum(Dn)>=-1e-8*maximum(Dn)
		@test count(Dn.<1e-6*maximum(Dn))>0          #part of the crack is closed
		DnB[name]=Dn
	end

	#C: two cracks
	Shift=[4.0 0 0]
	P1c=[P1;P1.+Shift]; P2c=[P2;P2.+Shift]; P3c=[P3;P3.+Shift]
	(FaceNormalVectorC,MidPointC)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(P1c,P2c,P3c)
	ZeroC=zeros(2n)
	BCC=MixedBoundaryConditionsFluidVolume(MixedBoundaryConditions(Stresses(ZeroC,ZeroC,fill(-1.0,2n),ZeroC,ZeroC,ZeroC),Tractions(ZeroC,ZeroC,ZeroC)),[0.5*V1,1.5*V1])
	(DnC,~,~)=CutAndDisplaceJulia.SlipCalculator3D(P1c,P2c,P3c,ν,G,λ,MidPointC,FaceNormalVectorC,HSFlag,BCC,[fill(1,n);fill(2,n)])
	AreaC=[Area;Area]
	@test abs(sum(DnC[1:n].*Area)/(0.5*V1)-1)<1e-8
	@test abs(sum(DnC[n+1:2n].*Area)/(1.5*V1)-1)<1e-8

	#Solver level: default (reduced) and full system agree and satisfy the LCP
	for (seed,Tn0,grad) in ((1,-1.0,3.0),(2,-2.0,4.0),(3,0.5,-2.0),(4,-3.0,6.0))
		(a,b,g,w,V)=FluidTestLCP(150,seed,Tn0,grad,0.5)
		(x,yy,P,Its,Conv)=CutAndDisplaceJulia.FluidVolumeFischerNewton(a,b,g,w,V)
		(xF,yF,PF,ItsF,ConvF)=CutAndDisplaceJulia.FluidVolumeNewton(a,b,g,w,V,false)
		@test Conv && ConvF
		Tol=1e-8*max(maximum(abs.(b)),1.0)
		for (xs,ys) in ((x,yy),(xF,yF))
			@test minimum(xs)>=0
			@test minimum(ys)>=-Tol
			@test maximum(abs.(xs.*ys))<=Tol*max(maximum(xs),1.0)
			@test maximum(abs.(w'*ys.-V))<=1e-8*maximum(V)
		end
		#x itself is not unique (see the top of FluidPressureSolver.jl: on open elements x4/y2 and x5/y3 can move
		#together), so compare what is: opening, slip (y-x of rows 2,3 and 4,5) and contact tractions (x rows 1-3)
		m=150; L1=1:m; L2=m+1:2m; L3=2m+1:3m; L4=3m+1:4m; L5=4m+1:5m
		Phys(x,y)=[y[L1]; y[L2]-x[L4]; y[L3]-x[L5]; x[1:3m]]
		@test maximum(abs.(Phys(x,yy).-Phys(xF,yF)))<=1e-6*max(maximum(abs.(Phys(xF,yF))),1.0)
		@test maximum(abs.(P.-PF))<=1e-8*maximum(abs.(PF))
		@test count(x[1:150].>0)>0                    #part of the crack is in contact
	end

	if MakeFigures
		fig=CM.Figure(size=(1300,550))
		ax1=CM.Axis(fig[1,1],xlabel="opening, unit tension",ylabel="opening, fluid volume solve",aspect=CM.DataAspect(),
			title="Case A (fully open): pressure 2 against compression 1\nmax |difference| = $(round(maximum(abs.(DnA.-DnRef)),sigdigits=2))")
		CM.scatter!(ax1,DnRef,DnA,markersize=6,color=:blue)
		CM.lines!(ax1,[0,maximum(DnRef)],[0,maximum(DnRef)],color=:black)
		ax2=CM.Axis(fig[1,2],xlabel="y (compression increases with -y)",ylabel="opening Dn",
			title="Cases B (0.3 V1) and B2 (0.6 V1): partly closed crack")
		CM.scatter!(ax2,y,DnB["B"],markersize=6,color=:red,label="B, 0.3 V1")
		CM.scatter!(ax2,y,DnB["B2"],markersize=6,color=:blue,label="B2, 0.6 V1")
		CM.axislegend(ax2,position=:lt)
		CM.Label(fig[0,:],"Unit test 12: fluid volume solve",fontsize=20)
		SaveTestFigure(fig,"ut_12_FluidSolverReduced")
	end
end
