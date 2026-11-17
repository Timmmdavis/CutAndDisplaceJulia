# Unit test 6: fluid-filled cracks with a fixed volume, SlipCalculator3D(...,MixedBoundaryConditionsFluidVolume,...).
# The fluid pressure is solved for together with the friction/contact problem (FluidPressureSolver.jl).
#  A) Penny under uniform compression σzz = -1, asked to hold the volume V1 that a pressure of +1
#     (no compression) gives. The solve must find a fluid pressure of 2 (1 to cancel the
#     compression + 1), so the opening must equal the plain "pressure 1" solution exactly.
#  B) Compression that increases downwards (σzz = -1 + 2y) and a small volume (0.3·V1): the lower
#     part of the crack must close (contact), with no interpenetration, and the volume is still held.
#  C) Two separate cracks in one solve, each holding its own volume.
# Figure: A against the plain solve and Sneddon/Eshelby; B opening map with closed elements; B
# opening against y; C the two cracks.
include(joinpath(@__DIR__,"UnitTestHelpers.jl"))

@testset "Fluid volume solve" begin

	(Points,Triangles,P1,P2,P3,FaceNormalVector,MidPoint)=LoadSmallPennyMesh()
	n=size(P1,1)
	(Area,~)=TriAreas(P1,P2,P3)
	G=ShearModulus(1.0); ν=PoissonsRatio(0.25)
	(K,E,λ,ν,G)=CutAndDisplaceJulia.ElasticConstantsCheck(G,ν)
	HSFlag=0
	Zero=zeros(n)
	FluidBC(σzz,Volume)=MixedBoundaryConditionsFluidVolume(
		MixedBoundaryConditions(Stresses(Zero,Zero,σzz,Zero,Zero,Zero),Tractions(Zero,Zero,Zero)),Volume)

	#Reference: pressure of 1, no fluid volume constraint
	(DnRef,~,~)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,Tractions(ones(n),Zero,Zero),zeros(n))
	V1=sum(DnRef.*Area)

	#A) uniform compression of 1, hold V1
	(DnA,DssA,DdsA)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,FluidBC(fill(-1.0,n),V1),ones(n))
	@test sum(DnA.*Area) ≈ V1 rtol=1e-8
	@test maximum(abs.(DnA.-DnRef))<1e-6*maximum(DnRef)
	@test maximum(abs.([DssA;DdsA]))<1e-8*maximum(DnRef)

	#B) compression increasing downwards, small volume: bottom closes
	y=MidPoint[:,2]
	VB=0.3*V1
	(DnB,DssB,DdsB)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,FluidBC(-1.0.+2.0.*y,VB),ones(n))
	ClosedB=DnB.<1e-3*maximum(DnB)
	@test sum(DnB.*Area) ≈ VB rtol=1e-8
	@test minimum(DnB)>-1e-6*maximum(DnB)                  #no interpenetration
	@test sum(ClosedB)>n/10                                  #a real closed region...
	@test mean(y[ClosedB])<-0.3                              #...at the bottom
	@test all(DnB[y.>0.5].>1e-3*maximum(DnB))                #top stays open

	#C) two cracks: the penny and a copy shifted 4 along x, different volumes
	Shift=[4.0 0 0]
	P1c=[P1;P1.+Shift]; P2c=[P2;P2.+Shift]; P3c=[P3;P3.+Shift]
	(FaceNormalVectorC,MidPointC)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(P1c,P2c,P3c)
	Flag=[fill(1,n);fill(2,n)]
	VC=[0.5*V1,1.5*V1]
	ZeroC=zeros(2n)
	BCC=MixedBoundaryConditionsFluidVolume(MixedBoundaryConditions(Stresses(ZeroC,ZeroC,fill(-1.0,2n),ZeroC,ZeroC,ZeroC),Tractions(ZeroC,ZeroC,ZeroC)),VC)
	(DnC,~,~)=CutAndDisplaceJulia.SlipCalculator3D(P1c,P2c,P3c,ν,G,λ,MidPointC,FaceNormalVectorC,HSFlag,BCC,Flag)
	AreaC=[Area;Area]
	@test sum(DnC[Flag.==1].*AreaC[Flag.==1]) ≈ VC[1] rtol=1e-8
	@test sum(DnC[Flag.==2].*AreaC[Flag.==2]) ≈ VC[2] rtol=1e-8

	if MakeFigures
		r=hypot.(MidPoint[:,1],MidPoint[:,2])
		rr=range(0,maximum(r),length=200)
		Rim=maximum(hypot.(Points[:,2],Points[:,3]))
		un=CutAndDisplaceJulia.Eshelby1957_PennyCrackSlipProfile(G,ν,1.0,0.0,collect(rr),Rim)
		fig=CM.Figure(size=(1900,600))
		ax1=CM.Axis(fig[1,1],xlabel="distance from centre r",ylabel="opening Dn",
			title="A) σzz = −1 everywhere, hold volume V1\nV1 = volume for a pressure of 1 (no σzz)\nsolved opening = pressure-1 opening (max diff $(round(maximum(abs.(DnA.-DnRef)),sigdigits=2)))")
		CM.lines!(ax1,rr,un,color=:black,label="analytical (Eshelby/Sneddon), pressure 1")
		CM.scatter!(ax1,r,DnRef,color=(:gray,0.6),markersize=10,label="SlipCalculator3D, pressure 1")
		CM.scatter!(ax1,r,DnA,color=:blue,markersize=4,label="fluid volume solve")
		CM.axislegend(ax1,position=:lb,labelsize=11)
		ax2=CM.Axis(fig[1,2],aspect=CM.DataAspect(),xlabel="x",ylabel="y",
			title="B) σzz = −1 + 2y, volume 0.3·V1\nopening; black outline: closed (in contact)\nvolume held: $(round(sum(DnB.*Area),sigdigits=6)) / $(round(VB,sigdigits=6))")
		m=FillTris!(ax2,P1,P2,P3,DnB,colormap=:viridis,colorrange=(0,maximum(DnB)))
		DrawTriEdges!(ax2,P1[ClosedB,:],P2[ClosedB,:],P3[ClosedB,:],color=:black,linewidth=0.8)
		CM.Colorbar(fig[1,3],m,label="opening Dn")
		ax3=CM.Axis(fig[1,4],xlabel="y (up the crack)",ylabel="opening Dn",
			title="B) opening against y\n$(sum(ClosedB)) of $n elements closed\nmin Dn = $(round(minimum(DnB),sigdigits=2)) (no interpenetration)\n(spread at one y: elements at different x)")
		CM.scatter!(ax3,y[.!ClosedB],DnB[.!ClosedB],color=:blue,markersize=5,label="open")
		CM.scatter!(ax3,y[ClosedB],DnB[ClosedB],color=:red,markersize=5,label="closed")
		CM.hlines!(ax3,[0],color=:black)
		CM.axislegend(ax3,position=:lt)
		ax4=CM.Axis(fig[1,5],aspect=CM.DataAspect(),xlabel="x",ylabel="y",
			title="C) two cracks in one solve, σzz = −1\nleft: $(round(sum(DnC[Flag.==1].*AreaC[Flag.==1])/V1,digits=6))·V1 (asked 0.5)\nright: $(round(sum(DnC[Flag.==2].*AreaC[Flag.==2])/V1,digits=6))·V1 (asked 1.5)")
		m4=FillTris!(ax4,P1c,P2c,P3c,DnC,colormap=:viridis,colorrange=(0,maximum(DnC)))
		CM.Colorbar(fig[1,6],m4,label="opening Dn",height=CM.Relative(0.6))
		CM.colsize!(fig.layout,5,CM.Auto(2))
		CM.Label(fig[0,:],"Unit test 6: fluid-filled cracks holding a fixed volume (pressure solved with the contact problem)",fontsize=20)
		SaveTestFigure(fig,"ut_06_FluidVolumeSolve")
	end
end
