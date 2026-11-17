# Unit test 4: remeshing with CGAL (needs BuildCGAL; skipped inside Pkg.test).
# The three CGAL tools the propagation loop uses, on the penny mesh:
#  - PolygonRemeshingCGAL: isotropic remesh to a target edge length (from GetDesiredEdgeLength for
#    300 triangles). Checks: typical edge length near the target, about 300 triangles, same area,
#    rim not moved outwards.
#  - AdvancingFrontCGAL: rebuilds a surface from points only (the crack points after a propagation
#    step). Checks: the rebuilt surface covers the same area.
#  - ConnectedComponentsCGAL + ConnectedComponentsReader: label separate pieces of a mesh. Checks:
#    the penny plus a shifted copy gives exactly two labels, one per piece.
# Figure: input, remeshed, edge-length histogram vs target, rebuilt-from-points, two components.
include(joinpath(@__DIR__,"UnitTestHelpers.jl"))

#Lengths of all triangle edges (each interior edge appears twice, fine for a histogram)
EdgeLengths(P1,P2,P3)=vec([sqrt.(sum((P1.-P2).^2,dims=2)) sqrt.(sum((P2.-P3).^2,dims=2)) sqrt.(sum((P3.-P1).^2,dims=2))])

if !isdefined(Main,:BuildCGAL)
	@info "BuildCGAL not loaded: skipping the CGAL remeshing test (run RunUnitTests.jl to include it)"
else
@testset "Remeshing (CGAL)" begin

	(Points,Triangles,P1,P2,P3,FaceNormalVector,MidPoint)=LoadPennyMesh()
	(Area0,~)=TriAreas(P1,P2,P3)
	R0=maximum(hypot.(Points[:,2],Points[:,3]))
	Dir=mktempdir()

	#Isotropic remeshing to the size that gives ~300 triangles
	NoTris=300
	(target_edge_length,~)=CutAndDisplaceJulia.GetDesiredEdgeLength(P1,P2,P3,NoTris)
	(Pts1,Tris1)=cd(Dir) do
		f=CutAndDisplaceJulia.OFFExport(Points,Triangles,size(Triangles,1),size(Points,1),"BeforeRemeshing")
		CutAndDisplaceJulia.OFFReader(BuildCGAL.PolygonRemeshingCGAL(f,target_edge_length))
	end
	(Q1,Q2,Q3)=CutAndDisplaceJulia.CreateP1P2P3(Tris1,Pts1)
	(Area1,~)=TriAreas(Q1,Q2,Q3)
	L=EdgeLengths(Q1,Q2,Q3)
	@test median(L) ≈ target_edge_length rtol=0.25
	@test 0.6*NoTris < size(Tris1,1) < 1.6*NoTris
	@test sum(Area1) ≈ sum(Area0) rtol=0.01
	@test maximum(hypot.(Pts1[:,2],Pts1[:,3]))<=R0+1e-6

	#Advancing front: rebuild a surface from the remeshed points alone
	(Pts2,Tris2)=cd(Dir) do
		f=CutAndDisplaceJulia.xyzExport(Pts1[:,2],Pts1[:,3],Pts1[:,4],"PointsOnly")
		CutAndDisplaceJulia.OFFReader(BuildCGAL.AdvancingFrontCGAL(f))
	end
	(S1,S2,S3)=CutAndDisplaceJulia.CreateP1P2P3(Tris2,Pts2)
	(Area2,~)=TriAreas(S1,S2,S3)
	@test sum(Area2) ≈ sum(Area1) rtol=0.03

	#Connected components: the remeshed penny and a copy shifted 3 along x, in one file
	Shift=3.0
	PtsTwo=[Pts1; [Pts1[:,1].+size(Pts1,1) Pts1[:,2].+Shift Pts1[:,3] Pts1[:,4]]]
	TrisTwo=[Tris1; Tris1.+size(Pts1,1)]
	Flags=cd(Dir) do
		f=CutAndDisplaceJulia.OFFExport(PtsTwo,TrisTwo,size(TrisTwo,1),size(PtsTwo,1),"TwoPieces")
		CutAndDisplaceJulia.ConnectedComponentsReader(BuildCGAL.ConnectedComponentsCGAL(f))
	end
	n1=size(Tris1,1)
	@test length(Flags)==size(TrisTwo,1)
	@test sort(unique(Flags))==[1,2]
	@test length(unique(Flags[1:n1]))==1 && length(unique(Flags[n1+1:end]))==1
	@test Flags[1]!=Flags[end]

	if MakeFigures
		fig=CM.Figure(size=(2000,520))
		ax1=CM.Axis(fig[1,1],aspect=CM.DataAspect(),xlabel="x",ylabel="y",
			title="Input penny\n$(size(Triangles,1)) triangles, area $(round(sum(Area0),digits=4))")
		DrawTriEdges!(ax1,P1,P2,P3)
		ax2=CM.Axis(fig[1,2],aspect=CM.DataAspect(),xlabel="x",ylabel="y",
			title="PolygonRemeshingCGAL\ntarget edge $(round(target_edge_length,digits=3)) (for $NoTris tris)\n→ $(size(Tris1,1)) triangles, area $(round(sum(Area1),digits=4))")
		DrawTriEdges!(ax2,Q1,Q2,Q3,color=(:black,0.6))
		ax3=CM.Axis(fig[1,3],xlabel="edge length",ylabel="number of edges",
			title="Remeshed edge lengths\nmedian $(round(median(L),digits=3)) vs target $(round(target_edge_length,digits=3))")
		CM.hist!(ax3,L,bins=30,color=(:gray,0.7))
		CM.vlines!(ax3,[target_edge_length],color=:red,linewidth=2,label="target")
		CM.axislegend(ax3)
		ax4=CM.Axis(fig[1,4],aspect=CM.DataAspect(),xlabel="x",ylabel="y",
			title="AdvancingFrontCGAL from the points only\n$(size(Tris2,1)) triangles, area $(round(sum(Area2),digits=4))")
		DrawTriEdges!(ax4,S1,S2,S3,color=(:blue,0.6))
		CM.scatter!(ax4,Pts1[:,2],Pts1[:,3],color=:black,markersize=3)
		ax5=CM.Axis(fig[1,5],aspect=CM.DataAspect(),xlabel="x",ylabel="y",
			title="ConnectedComponentsCGAL\ntwo separate pennies → labels $(sort(unique(Flags)))\n(orange = label 1, teal = label 2)")
		(T1,T2,T3)=CutAndDisplaceJulia.CreateP1P2P3(TrisTwo,PtsTwo)
		m=FillTris!(ax5,T1,T2,T3,Float64.(Flags),colormap=CM.cgrad([:orange,:teal],2,categorical=true),colorrange=(0.5,2.5))
		DrawTriEdges!(ax5,T1,T2,T3,color=(:white,0.6),linewidth=0.4)
		CM.colsize!(fig.layout,5,CM.Auto(1.8))
		CM.Label(fig[0,:],"Unit test 4: CGAL remeshing tools used in the propagation loop",fontsize=20)
		SaveTestFigure(fig,"ut_04_RemeshingCGAL")
	end
end
end
