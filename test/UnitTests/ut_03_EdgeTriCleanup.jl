# Unit test 3: cleaning up the triangles on the crack rim (tip triangles).
# Remeshing can leave thin "sliver" triangles on the rim, which give poor stress intensities.
# The propagation loop fixes these with two functions, tested here on the penny mesh (which has
# some rim slivers, see unit test 2):
#  - CleanEdgeTris: removes slivers that hang off the rim, and merges rim triangles that share
#    one inner point. Checks: the worst rim angle improves a lot and little area is lost.
#  - IsosceliseEdgeTris: moves the inner point of each rim triangle so the triangle becomes
#    closer to isosceles. Checks: the area is unchanged, the rim does not move, and rim triangles
#    get fatter (larger smallest-angle).
# Figure: the mesh at each stage coloured by smallest angle, rim triangles outlined, and the
# distribution of rim-triangle angles at each stage.
include(joinpath(@__DIR__,"UnitTestHelpers.jl"))

#Rim (tip) triangles and the smallest internal angle of every triangle
function RimTrianglesAndAngles(P1,P2,P3,FaceNormalVector,MidPoint)
	(F1,F2,F3)=CutAndDisplaceJulia.GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector)
	Rim=F1.FreeFlg .| F2.FreeFlg .| F3.FreeFlg
	(a,b,c)=CutAndDisplaceJulia.CalculateInternalTriAngles(P1,P2,P3)
	return Rim,min.(a,b,c)
end

@testset "Edge (rim) triangle cleanup" begin

	(Points,Triangles,P1,P2,P3,FaceNormalVector,MidPoint)=LoadPennyMesh()
	(Area0,~)=TriAreas(P1,P2,P3)
	(Rim0,MinAng0)=RimTrianglesAndAngles(P1,P2,P3,FaceNormalVector,MidPoint)
	R0=maximum(hypot.(P1[:,1],P1[:,2]))

	#Stage 1: remove rim slivers and merge rim triangles sharing an inner point
	(Q1,Q2,Q3,Tris1,Pts1,MidPoint1,FaceNormalVector1)=CutAndDisplaceJulia.CleanEdgeTris(Points,Triangles)
	(Area1,~)=TriAreas(Q1,Q2,Q3)
	(Rim1,MinAng1)=RimTrianglesAndAngles(Q1,Q2,Q3,FaceNormalVector1,MidPoint1)
	@test minimum(MinAng0[Rim0])<5                 #the input really has slivers
	@test minimum(MinAng1[Rim1])>10                #...and they are gone
	@test abs(sum(Area1)-sum(Area0))/sum(Area0)<0.01
	@test maximum(hypot.(Q1[:,1],Q1[:,2]))<=R0+1e-12

	#Stage 2: make the rim triangles closer to isosceles
	Before=deepcopy((Q1,Q2,Q3))
	(S1,S2,S3,Tris2,Pts2,MidPoint2,FaceNormalVector2)=CutAndDisplaceJulia.IsosceliseEdgeTris(MidPoint1,Q1,Q2,Q3,Tris1,Pts1,FaceNormalVector1)
	#all its outputs describe the same triangles (the propagation loop uses both Triangles and MidPoint)
	@test size(Tris2,1)==size(S1,1)==size(MidPoint2,1)==size(FaceNormalVector2,1)
	(Area2,~)=TriAreas(S1,S2,S3)
	(Rim2,MinAng2)=RimTrianglesAndAngles(S1,S2,S3,FaceNormalVector2,MidPoint2)
	@test sum(Area2) ≈ sum(Area1) rtol=1e-10        #inner points move, the outline does not
	@test (Q1,Q2,Q3)==Before                        #the input mesh is left unchanged
	#After each step the outline is the rim only: no boundary edge cuts into the mesh
	Notch(A,B,C)=[(a,b) for (a,b) in BoundaryEdgesBruteForce(A,B,C) if min(hypot(a[1],a[2]),hypot(b[1],b[2]))<0.9]
	@test isempty(Notch(Q1,Q2,Q3))
	@test isempty(Notch(S1,S2,S3))
	@test length(BoundaryEdgesBruteForce(S1,S2,S3))==length(BoundaryEdgesBruteForce(Q1,Q2,Q3))
	@test mean(MinAng2[Rim2])>mean(MinAng1[Rim1])+3
	@test sum(MinAng2[Rim2].<20)<sum(MinAng1[Rim1].<20)
	@test all(Area2.>0)

	if MakeFigures
		fig=CM.Figure(size=(1700,620))
		Stages=(("Input mesh",P1,P2,P3,Rim0,MinAng0),
		        ("After CleanEdgeTris",Q1,Q2,Q3,Rim1,MinAng1),
		        ("After IsosceliseEdgeTris",S1,S2,S3,Rim2,MinAng2))
		local m
		for (k,(name,A,B,C,Rim,MinAng)) in enumerate(Stages)
			ax=CM.Axis(fig[1,k],aspect=CM.DataAspect(),xlabel="x",ylabel="y",
				title="$name\n$(size(A,1)) triangles, $(sum(Rim)) on the rim\nworst rim angle $(round(minimum(MinAng[Rim]),digits=1))°, mean $(round(mean(MinAng[Rim]),digits=1))°")
			m=FillTris!(ax,A,B,C,MinAng,colormap=:viridis,colorrange=(0,60))
			DrawTriEdges!(ax,A,B,C,color=(:white,0.5),linewidth=0.4)
			DrawTriEdges!(ax,A[Rim,:],B[Rim,:],C[Rim,:],color=:red,linewidth=1.0)
		end
		CM.Colorbar(fig[1,4],m,label="smallest internal angle [°]")
		ax=CM.Axis(fig[1,5],xlabel="smallest angle of rim triangle [°]",ylabel="fraction of rim triangles",
			title="Rim triangles get fatter\n(higher is better, 60° = equilateral)")
		for (k,(name,A,B,C,Rim,MinAng)) in enumerate(Stages)
			v=sort(MinAng[Rim])
			CM.stairs!(ax,v,(1:length(v))./length(v),label=name,linewidth=2,color=(:gray,:orange,:blue)[k])
		end
		CM.axislegend(ax,position=:rb,labelsize=11)
		CM.Label(fig[2,1:3],"red outlines: rim (tip) triangles",tellwidth=false)
		CM.Label(fig[0,:],"Unit test 3: cleaning up rim triangles before computing stress intensity",fontsize=20)
		SaveTestFigure(fig,"ut_03_EdgeTriCleanup")
	end
end
