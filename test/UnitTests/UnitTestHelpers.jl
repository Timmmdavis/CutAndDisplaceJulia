# Shared set-up for the unit tests in test/UnitTests.
#
# Each unit test checks one part of the code on a small problem (seconds, not minutes), using
# @test so it can run inside Pkg.test (runtests.jl). When run through RunUnitTests.jl it also
# saves a labelled figure to test/UnitTests/figures/<test name>.png so you can SEE what was tested.
#
# Figures need CairoMakie. RunUnitTests.jl loads it (from CutAndDisplaceJuliaPlots) and binds it
# to the name CM before including the tests; inside Pkg.test, CM is not defined and no figures
# are drawn.

#Only define these once, however many unit tests include this file
if !@isdefined(UnitTestHelpersLoaded)
UnitTestHelpersLoaded=true

using Test
using Statistics
using LinearAlgebra
using CutAndDisplaceJulia

MakeFigures = isdefined(Main, :CM)
FigureDir = joinpath(@__DIR__, "figures")
MakeFigures && mkpath(FigureDir)

#Unit radius penny-shaped crack in the x-y plane (z = 0), ~700 triangles
function LoadPennyMesh()
	SurfaceDir=CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"CircleMesh_1a_500Faces.off")
	(Points,Triangles)=CutAndDisplaceJulia.OFFReader(SurfaceDir)
	(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3(Triangles,Points)
	(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
	return Points,Triangles,P1,P2,P3,FaceNormalVector,MidPoint
end

#Smaller unit-radius penny (~270 near-equilateral triangles) in the x-y plane, for the solver tests
#(the stored .stl is vertical: swap its x and z so the normals point along z)
function LoadSmallPennyMesh()
	SurfaceDir=CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"VertPenny-300-EqEdges.stl")
	(Points,Triangles)=CutAndDisplaceJulia.STLReader(SurfaceDir)
	Points[:,[2,4]]=Points[:,[4,2]]
	Points[:,4].=0.0
	(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3(Triangles,Points)
	(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
	return Points,Triangles,P1,P2,P3,FaceNormalVector,MidPoint
end

#Triangle areas as a vector
function TriAreas(P1,P2,P3)
	(Area,HalfPerimeter)=CutAndDisplaceJulia.AreaOfTriangle3D(P1[:,1],P1[:,2],P1[:,3],P2[:,1],P2[:,2],P2[:,3],P3[:,1],P3[:,2],P3[:,3])
	return vec(Area),vec(HalfPerimeter)
end

#Edges used by only one triangle (the mesh boundary), found independently of the package code.
#Returns a list of (point a, point b) as 1×3 rows.
function BoundaryEdgesBruteForce(P1,P2,P3; digits=10)
	key(p)=Tuple(round.(p,digits=digits))
	count=Dict{Any,Int}()
	ends=Dict{Any,Any}()
	for i=1:size(P1,1)
		for (a,b) in ((P1[i,:],P2[i,:]),(P2[i,:],P3[i,:]),(P3[i,:],P1[i,:]))
			k=Set([key(a),key(b)])
			count[k]=get(count,k,0)+1
			ends[k]=(a,b)
		end
	end
	return [ends[k] for k in keys(count) if count[k]==1]
end

#Save a figure (only when figures are on) and say where it went
function SaveTestFigure(fig,name)
	MakeFigures || return nothing
	f=joinpath(FigureDir,"$name.png")
	CM.save(f,fig)
	println("    figure: $f")
	return f
end

#Draw triangle outlines on a 2D axis (map view by default: x across, y up)
function DrawTriEdges!(ax,P1,P2,P3; ia=1, ib=2, color=(:black,0.35), linewidth=0.6)
	seg=CM.Point2f[]
	for i=1:size(P1,1), (a,b) in ((P1,P2),(P2,P3),(P3,P1))
		push!(seg,CM.Point2f(a[i,ia],a[i,ib])); push!(seg,CM.Point2f(b[i,ia],b[i,ib]))
	end
	CM.linesegments!(ax,seg; color=color, linewidth=linewidth)
end

#Fill triangles with one colour value each (map view by default)
function FillTris!(ax,P1,P2,P3,Values; ia=1, ib=2, colormap=:viridis, colorrange=nothing)
	n=size(P1,1)
	V=Matrix{Float64}(undef,3n,2)
	for i=1:n
		V[3i-2,:]=[P1[i,ia],P1[i,ib]]; V[3i-1,:]=[P2[i,ia],P2[i,ib]]; V[3i,:]=[P3[i,ia],P3[i,ib]]
	end
	F=reshape(collect(1:3n),3,n)'
	c=repeat(vec(Values),inner=3)
	cr=colorrange===nothing ? extrema(filter(isfinite,c)) : colorrange
	return CM.mesh!(ax,V,F; color=c, colormap=colormap, colorrange=cr, shading=CM.NoShading)
end

end
