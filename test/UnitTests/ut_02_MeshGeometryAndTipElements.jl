# Unit test 2: triangle geometry and finding the crack-tip (edge) elements.
# On the unit-radius penny mesh:
#  - AreaOfTriangle3D: areas add up to the area inside the rim (shoelace formula on the rim edges)
#  - CreateFaceNormalAndMidPoint: normals are unit vectors (±z here); "MidPoint" is the incentre
#    (corners weighted by the opposite edge lengths), checked against that formula
#  - CalculateInternalTriAngles: the three angles of every triangle add to 180°
#  - GetCrackTipElements3D: flags exactly the edges that belong to only one triangle (checked by
#    brute force), their midpoints sit on the rim, and the midpoint-to-edge vector points outwards
#  - GetDesiredEdgeLength: equilateral triangles of that size and number give the same total area
# Figure: mesh coloured by smallest internal angle; tip edges and their outward vectors; angle histogram.
include(joinpath(@__DIR__,"UnitTestHelpers.jl"))

@testset "Mesh geometry and tip elements" begin

	(Points,Triangles,P1,P2,P3,FaceNormalVector,MidPoint)=LoadPennyMesh()
	n=size(P1,1)

	#Areas
	(Area,HalfPerimeter)=TriAreas(P1,P2,P3)
	@test all(Area.>0)
	RimEdges=BoundaryEdgesBruteForce(P1,P2,P3)
	ShoelaceArea=abs(sum(a[1]*b[2]-b[1]*a[2] for (a,b) in RimEdges))/2   #rim edges are not ordered,
	@test sum(Area) ≈ ShoelaceArea rtol=1e-8                                  #but each adds its own wedge
	@test sum(Area) ≈ π rtol=0.05                                            #roughly a unit circle

	#Normals and midpoints
	@test all(abs.(sqrt.(sum(FaceNormalVector.^2,dims=2)).-1).<1e-10)
	@test all(abs.(abs.(FaceNormalVector[:,3]).-1).<1e-10)
	La=sqrt.(sum((P2.-P3).^2,dims=2)); Lb=sqrt.(sum((P1.-P3).^2,dims=2)); Lc=sqrt.(sum((P1.-P2).^2,dims=2))
	Incentre=(P1.*La.+P2.*Lb.+P3.*Lc)./(La.+Lb.+Lc)
	@test maximum(abs.(MidPoint.-Incentre))<1e-12

	#Internal angles (degrees)
	(AngA,AngB,AngC)=CutAndDisplaceJulia.CalculateInternalTriAngles(P1,P2,P3)
	@test maximum(abs.(AngA.+AngB.+AngC.-180))<1e-8
	MinAngle=min.(AngA,AngB,AngC)          #(this mesh has a few slivers on the rim: see unit test 4)

	#Tip elements: each triangle has 3 connections (P1P2, P1P3, P2P3); FreeFlg marks edges on the rim
	(FeP1P2S,FeP1P3S,FeP2P3S)=CutAndDisplaceJulia.GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector)
	Fe=(FeP1P2S,FeP1P3S,FeP2P3S)
	NoTipEdges=sum(sum(F.FreeFlg) for F in Fe)
	@test NoTipEdges==length(RimEdges)
	TipMid=vcat([F.FeMd[F.FreeFlg,:] for F in Fe]...)
	TipOut=vcat([F.FeM2Ev[F.FreeFlg,:] for F in Fe]...)          #triangle midpoint -> edge
	TipTri=vcat([MidPoint[F.FreeFlg,:] for F in Fe]...)          #midpoint of the triangle it belongs to
	TipLen=vcat([F.FeLe[F.FreeFlg] for F in Fe]...)
	RimR=hypot.(TipMid[:,1],TipMid[:,2])
	@test all(RimR.>0.95) && all(RimR.<1.01)                    #edge midpoints on the rim
	@test all(sum(TipOut[:,1:2].*TipMid[:,1:2],dims=2).>0)      #pointing outwards
	@test all(abs.(sqrt.(sum(TipOut.^2,dims=2)).-1).<1e-10)     #unit vectors
	@test all(sum(TipOut.*(TipMid.-TipTri),dims=2).>0)          #from the triangle towards its edge
	@test sum(TipLen) ≈ 2π rtol=0.02                            #rim length ≈ circumference

	#Desired edge length: NoTris equilateral triangles with this edge cover the same area
	NoTris=300
	(target_edge_length,max_target_edge_length)=CutAndDisplaceJulia.GetDesiredEdgeLength(P1,P2,P3,NoTris)
	@test NoTris*(sqrt(3)/4)*target_edge_length^2 ≈ sum(Area) rtol=1e-10
	@test max_target_edge_length ≈ maximum(HalfPerimeter)*(2/3)

	if MakeFigures
		fig=CM.Figure(size=(1500,560))
		ax1=CM.Axis(fig[1,1],aspect=CM.DataAspect(),xlabel="x",ylabel="y",
			title="Penny mesh, $n triangles, coloured by smallest angle\nΣ triangle areas = $(round(sum(Area),digits=6))\narea inside rim (shoelace) = $(round(ShoelaceArea,digits=6))")
		m=FillTris!(ax1,P1,P2,P3,MinAngle,colormap=:viridis,colorrange=(20,60))
		DrawTriEdges!(ax1,P1,P2,P3,color=(:white,0.5),linewidth=0.4)
		CM.Colorbar(fig[1,2],m,label="smallest internal angle [°]")
		ax2=CM.Axis(fig[1,3],aspect=CM.DataAspect(),xlabel="x",ylabel="y",
			title="GetCrackTipElements3D (zoom on one quadrant)\n$NoTipEdges tip edges, brute force count $(length(RimEdges))\narrows: tip triangle midpoint → edge midpoint",
			limits=(0.3,1.05,0.3,1.05))
		DrawTriEdges!(ax2,P1,P2,P3,color=(:black,0.15))
		for (a,b) in RimEdges
			CM.lines!(ax2,[a[1],b[1]],[a[2],b[2]],color=:red,linewidth=3)
		end
		CM.scatter!(ax2,TipTri[:,1],TipTri[:,2],color=:blue,markersize=5)
		CM.arrows2d!(ax2,TipTri[:,1],TipTri[:,2],TipMid[:,1].-TipTri[:,1],TipMid[:,2].-TipTri[:,2],color=:blue,shaftwidth=1.5,tipwidth=8,tiplength=8)
		CM.lines!(ax2,[NaN],[NaN],color=:red,linewidth=3,label="edges used by one triangle")
		CM.lines!(ax2,[NaN],[NaN],color=:blue,label="tip triangle → its free edge")
		CM.axislegend(ax2,position=:rb,framevisible=true)
		ax3=CM.Axis(fig[1,4],xlabel="internal angle [°]",ylabel="number of angles",
			title="CalculateInternalTriAngles (all $(3n) angles)\nnear 0° and 180°: slivers on the rim\nmax |A+B+C−180°| = $(round(maximum(abs.(AngA.+AngB.+AngC.-180)),sigdigits=2))°")
		CM.hist!(ax3,[AngA;AngB;AngC],bins=40,color=(:gray,0.7))
		CM.vlines!(ax3,[60],color=:black,linestyle=:dash)
		CM.text!(ax3,61,0,text="60° (equilateral)",align=(:left,:bottom))
		CM.Label(fig[0,:],"Unit test 2: triangle geometry and crack-tip elements on a unit penny",fontsize=20)
		SaveTestFigure(fig,"ut_02_MeshGeometryAndTipElements")
	end
end
