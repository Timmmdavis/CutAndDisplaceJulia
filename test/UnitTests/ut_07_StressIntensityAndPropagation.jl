# Unit test 7: stress intensity at the crack tip and one propagation step.
# Small penny (unit radius) opened by a pressure that increases upwards, Tn = 1 + 0.5y.
#  - StressIntensity3D: KI around the rim matches the analytical penny solution
#      KI = 2·p0·sqrt(a/π) + (4/(3π))·g·a·sqrt(πa)·(y/a)     (uniform part + gradient part, Tada)
#  - MovingAverageOfStressIntensity: smooths K along the rim (smaller jumps between neighbours,
#    same mean, no new extremes)
#  - PropagateFracture: only edges with K above KCrit grow; each new point lies outside the rim,
#    in the crack plane (pure opening), at most one element length away, with the largest K
#    moving the full length
#  - With BuildCGAL: AdvancingFrontCGAL + RemoveDodgyNewEdges rebuild the larger crack; it grows
#    at the top (high K) and not at the bottom
# Figure: KI around the rim against the analytical curve; the step (old rim, tip K/Kc, new points,
# and the rebuilt mesh when CGAL is available).
include(joinpath(@__DIR__,"UnitTestHelpers.jl"))

@testset "Stress intensity and one propagation step" begin

	(Points,Triangles,P1,P2,P3,FaceNormalVector,MidPoint)=LoadSmallPennyMesh()
	n=size(P1,1)
	(Area,HalfPerimeter)=TriAreas(P1,P2,P3)
	a=maximum(hypot.(Points[:,2],Points[:,3]))
	G=ShearModulus(1.0); ν=PoissonsRatio(0.25)
	(K,E,λ,ν,G)=CutAndDisplaceJulia.ElasticConstantsCheck(G,ν)

	#Opening under Tn = p0 + g*y
	p0=1.0; g=0.5
	Tn=p0.+g.*MidPoint[:,2]
	(Dn,Dss,Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,0,Tractions(Tn,zeros(n),zeros(n)),zeros(n))

	#Stress intensity on the tip elements
	(FeP1P2S,FeP1P3S,FeP2P3S)=CutAndDisplaceJulia.GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector)
	(FeP1P2S,FeP1P3S,FeP2P3S)=CutAndDisplaceJulia.StressIntensity3D(Dn,Dss,Dds,G,ν,FaceNormalVector,FeP1P2S,FeP1P3S,FeP2P3S)
	Fe=(FeP1P2S,FeP1P3S,FeP2P3S)
	TipXYZ=vcat([F.FeMd[F.FreeFlg,:] for F in Fe]...)
	KI=vcat([F.K1[F.FreeFlg] for F in Fe]...)
	KIAn=2*p0*sqrt(a/π).+(4/(3π))*g*a*sqrt(π*a).*(TipXYZ[:,2]./a)
	θ=atan.(TipXYZ[:,2],TipXYZ[:,1])
	@test all(isfinite.(KI))
	@test median(abs.(KI.-KIAn)./KIAn)<0.10                  #typically within a few %
	@test cor(KI,KIAn)>0.9                                    #follows the gradient

	#Smoothing along the rim
	Fe0=deepcopy(Fe)
	(S1,S2,S3)=CutAndDisplaceJulia.MovingAverageOfStressIntensity(3,P1,P2,P3,FaceNormalVector,MidPoint,deepcopy(FeP1P2S),deepcopy(FeP1P3S),deepcopy(FeP2P3S))
	KEq0=vcat([F.StrainEnergy[F.FreeFlg] for F in Fe0]...)
	KEqS=vcat([F.StrainEnergy[F.FreeFlg] for F in (S1,S2,S3)]...)
	Roughness(K)=std(diff(K[sortperm(θ)]))                   #jumps between neighbours round the rim
	@test Roughness(KEqS)<Roughness(KEq0)
	@test minimum(KEq0)-1e-12<=minimum(KEqS) && maximum(KEqS)<=maximum(KEq0)+1e-12
	@test mean(KEqS) ≈ mean(KEq0) rtol=0.05

	#One propagation step: grow where K > KCrit
	KCrit=median(KEqS)
	AvgTriangleEdgeLength=mean(HalfPerimeter)*(2/3)
	(PV1,PV2,PV3,Still12,Still13,Still23,~,~,~)=CutAndDisplaceJulia.PropagateFracture(S1,S2,S3,FaceNormalVector,G,ν,KCrit,AvgTriangleEdgeLength)
	NewPts=Float64.([PV1;PV2;PV3])
	Grows=[F.StrainEnergy[i]>KCrit for F in (S1,S2,S3) for i in findall(F.FreeFlg)]
	@test size(NewPts,1)==sum(Grows)
	@test sum(Still12)+sum(Still13)+sum(Still23)==length(Grows)-sum(Grows)
	RNew=hypot.(NewPts[:,1],NewPts[:,2])
	@test all(RNew.>minimum(hypot.(TipXYZ[Grows,1],TipXYZ[Grows,2])))   #outside the rim
	@test maximum(abs.(NewPts[:,3]))<1e-8                                  #stays in the crack plane
	@test mean(NewPts[:,2])>0.3                                            #growth is at the top
	#distance moved: K/Kmax × element length, from the tip edge midpoint
	Step=[minimum(hypot.(TipXYZ[Grows,1].-p[1],TipXYZ[Grows,2].-p[2])) for p in eachrow(NewPts)]
	@test maximum(Step) ≈ AvgTriangleEdgeLength rtol=0.05
	@test all(Step.<=AvgTriangleEdgeLength*1.0001)

	#Rebuild the larger crack from points (needs BuildCGAL), as the propagation loop does
	HaveCGAL=isdefined(Main,:BuildCGAL)
	if HaveCGAL
		Dir=mktempdir()
		#rim points that did not grow + all mesh points + new points
		AllPts=unique([NewPts; P1; P2; P3],dims=1)
		(PtsN,TrisN)=cd(Dir) do
			f=CutAndDisplaceJulia.xyzExport(AllPts[:,1],AllPts[:,2],AllPts[:,3],"NewFracturePoints")
			CutAndDisplaceJulia.OFFReader(BuildCGAL.AdvancingFrontCGAL(f))
		end
		(FNN,MPN)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(PtsN,TrisN)
		(Q1,Q2,Q3)=CutAndDisplaceJulia.CreateP1P2P3(TrisN,PtsN)
		max_target_edge_length=maximum(HalfPerimeter)*(2/3)
		(Q1,Q2,Q3,PtsN,TrisN)=CutAndDisplaceJulia.RemoveDodgyNewEdges(Q1,Q2,Q3,PtsN,TrisN,FNN,MPN,NewPts,max_target_edge_length)
		(AreaN,~)=TriAreas(Q1,Q2,Q3)
		OldY=[P1[:,2];P2[:,2];P3[:,2]]
		@test sum(AreaN)>sum(Area)                                        #it grew...
		@test maximum(PtsN[:,3])>maximum(OldY)+0.3*AvgTriangleEdgeLength  #...at the top
		@test minimum(PtsN[:,3])≈minimum(OldY) atol=0.3*AvgTriangleEdgeLength   #...not at the bottom
		@test maximum(abs.(PtsN[:,4]))<1e-8                              #still flat (z = 0)
	end

	if MakeFigures
		fig=CM.Figure(size=(1500,620))
		o=sortperm(θ)
		ax1=CM.Axis(fig[1,1],xlabel="angle round the rim θ [°] (90° = top, +y)",ylabel="KI",xticks=-180:90:180,
			title="StressIntensity3D, Tn = $p0 + $g·y\nline: analytical penny (uniform + gradient)\nmedian error $(round(100*median(abs.(KI.-KIAn)./KIAn),digits=1))%")
		CM.lines!(ax1,rad2deg.(θ[o]),KIAn[o],color=:black,linewidth=2,label="analytical")
		CM.scatter!(ax1,rad2deg.(θ),KI,color=:blue,markersize=7,label="KI per tip element")
		CM.scatter!(ax1,rad2deg.(θ),KEqS,color=:orange,marker=:xcross,markersize=8,label="after MovingAverage (K equivalent)")
		CM.hlines!(ax1,[KCrit],color=:red,linestyle=:dash,label="KCrit used for the step")
		CM.axislegend(ax1,position=:lt,labelsize=11)
		ax2=CM.Axis(fig[1,2],aspect=CM.DataAspect(),xlabel="x",ylabel="y",
			title="PropagateFracture with KCrit = median K\n$(sum(Grows)) of $(length(Grows)) tip edges grow (K > KCrit)\nlongest step = element length $(round(AvgTriangleEdgeLength,digits=3))")
		DrawTriEdges!(ax2,P1,P2,P3,color=(:black,0.2))
		KoverKc=KEqS./KCrit
		sc=CM.scatter!(ax2,TipXYZ[:,1],TipXYZ[:,2],color=KoverKc,colormap=:coolwarm,colorrange=(0.5,1.5),markersize=10)
		CM.Colorbar(fig[1,3],sc,label="tip K / KCrit (red > 1 grows)")
		CM.scatter!(ax2,NewPts[:,1],NewPts[:,2],color=:black,marker=:star5,markersize=10,label="new tip points")
		CM.axislegend(ax2,position=:lb,labelsize=11)
		if HaveCGAL
			ax3=CM.Axis(fig[1,4],aspect=CM.DataAspect(),xlabel="x",ylabel="y",
				title="Rebuilt crack\n(AdvancingFrontCGAL + RemoveDodgyNewEdges)\narea $(round(sum(Area),digits=3)) → $(round(sum(AreaN),digits=3))\ngrey: old crack")
			FillTris!(ax3,P1,P2,P3,zeros(n),colormap=[:lightgray,:lightgray],colorrange=(0,1))
			DrawTriEdges!(ax3,Q1,Q2,Q3,color=(:blue,0.7))
			CM.scatter!(ax3,NewPts[:,1],NewPts[:,2],color=:black,marker=:star5,markersize=8)
		end
		CM.Label(fig[0,:],"Unit test 7: stress intensity around the tip and one propagation step",fontsize=20)
		SaveTestFigure(fig,"ut_07_StressIntensityAndPropagation")
	end
end
