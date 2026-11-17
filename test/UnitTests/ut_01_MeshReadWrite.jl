# Unit test 1: reading and writing meshes.
# The same unit-radius circle is stored as .off and as GoCad .ts; a vertical penny as .stl.
# Checks: each reader gives a sensible mesh, and writing then reading back (.off, .stl, .xyz)
# returns the same geometry. Figure: the meshes as read, and the written-then-read copies on top.
include(joinpath(@__DIR__,"UnitTestHelpers.jl"))

@testset "Mesh read/write" begin

	Dir=mktempdir()

	#Read the same circle from two formats
	(PtsOff,TrisOff)=CutAndDisplaceJulia.OFFReader(CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"CircleMesh_1a_500Faces.off"))
	(PtsTs,TrisTs)=CutAndDisplaceJulia.GoCadAsciiReader(CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"CircleMesh_1a_500Faces.ts"))
	(PtsStl,TrisStl)=CutAndDisplaceJulia.STLReader(CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"VertPenny-300-EqEdges.stl"))

	for (Pts,Tris) in ((PtsOff,TrisOff),(PtsTs,TrisTs),(PtsStl,TrisStl))
		@test size(Pts,2)==4                            #index, x, y, z
		@test size(Tris,2)==3
		@test minimum(Tris)>=1 && maximum(Tris)<=size(Pts,1)
	end
	#Both circles are unit radius in the x-y plane
	@test maximum(hypot.(PtsOff[:,2],PtsOff[:,3])) ≈ 1.0 atol=0.01
	@test maximum(hypot.(PtsTs[:,2],PtsTs[:,3])) ≈ 1.0 atol=0.01
	@test maximum(abs.(PtsOff[:,4]))<1e-9

	#.off: write, read back, same points and triangles
	(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3(TrisOff,PtsOff)
	OffFile=cd(Dir) do
		CutAndDisplaceJulia.OFFExport(PtsOff,TrisOff,size(TrisOff,1),size(PtsOff,1),"RoundTrip")
	end
	@test isfile(OffFile)
	(Pts2,Tris2)=CutAndDisplaceJulia.OFFReader(OffFile)
	@test Tris2==TrisOff
	@test maximum(abs.(Pts2[:,2:4].-PtsOff[:,2:4]))<1e-12

	#.stl: write, read back, same triangle corners (STL stores each triangle's own corners)
	(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(PtsOff,TrisOff)
	StlFile=cd(Dir) do
		CutAndDisplaceJulia.STLExport(P1,P2,P3,FaceNormalVector,"RoundTrip")
	end
	@test isfile(StlFile)
	(Pts3,Tris3)=CutAndDisplaceJulia.STLReader(StlFile)
	(Q1,Q2,Q3)=CutAndDisplaceJulia.CreateP1P2P3(Tris3,Pts3)
	@test size(Q1,1)==size(P1,1)
	@test maximum(abs.([Q1 Q2 Q3].-[P1 P2 P3]))<1e-5         #STL is written in single-ish precision

	#.xyz: write a point list, read it back as numbers
	XyzFile=cd(Dir) do
		CutAndDisplaceJulia.xyzExport(PtsOff[:,2],PtsOff[:,3],PtsOff[:,4],"RoundTrip")
	end
	@test isfile(XyzFile)
	Xyz=CutAndDisplaceJulia.readdlm(XyzFile)
	@test size(Xyz,1)==size(PtsOff,1)
	@test maximum(abs.(Xyz[:,1:3].-PtsOff[:,2:4]))<1e-6

	if MakeFigures
		fig=CM.Figure(size=(1500,560))
		ax1=CM.Axis(fig[1,1],aspect=CM.DataAspect(),xlabel="x",ylabel="y",
			title="Circle read from .off ($(size(TrisOff,1)) triangles)\nred: written to .off and read back (lies on top)\nblue dots: points written to .xyz and read back")
		(A1,A2,A3)=CutAndDisplaceJulia.CreateP1P2P3(TrisOff,PtsOff)
		DrawTriEdges!(ax1,A1,A2,A3)
		(B1,B2,B3)=CutAndDisplaceJulia.CreateP1P2P3(Tris2,Pts2)
		DrawTriEdges!(ax1,B1,B2,B3,color=(:red,0.6),linewidth=0.4)
		CM.scatter!(ax1,Xyz[:,1],Xyz[:,2],color=:blue,markersize=4)
		ax2=CM.Axis(fig[1,2],aspect=CM.DataAspect(),xlabel="x",ylabel="y",
			title="Unit circle read from GoCad .ts\n(a different mesh: $(size(TrisTs,1)) triangles)")
		(C1,C2,C3)=CutAndDisplaceJulia.CreateP1P2P3(TrisTs,PtsTs)
		DrawTriEdges!(ax2,C1,C2,C3)
		ax3=CM.Axis(fig[1,3],aspect=CM.DataAspect(),xlabel="first coordinate",ylabel="third coordinate",
			title="Vertical penny read from .stl ($(size(TrisStl,1)) triangles)")
		(D1,D2,D3)=CutAndDisplaceJulia.CreateP1P2P3(TrisStl,PtsStl)
		#show the two coordinates with the largest spread (the plane of the penny)
		sp=[maximum(PtsStl[:,k])-minimum(PtsStl[:,k]) for k=2:4]
		ord=sortperm(sp,rev=true)
		DrawTriEdges!(ax3,D1,D2,D3,ia=ord[1],ib=ord[2])
		ax3.xlabel="xyz"[ord[1]:ord[1]]; ax3.ylabel="xyz"[ord[2]:ord[2]]
		CM.Label(fig[0,:],"Unit test 1: mesh readers and writers (.off, .ts, .stl, .xyz)",fontsize=20)
		SaveTestFigure(fig,"ut_01_MeshReadWrite")
	end
end
