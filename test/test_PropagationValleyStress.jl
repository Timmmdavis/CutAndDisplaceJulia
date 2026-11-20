# Crack loaded by the stresses below an incised valley (Savage et al. 1984): testPropagationValleyStress.
#
# Usage (after using CutAndDisplaceJulia, CutAndDisplaceJuliaPlots, BuildCGAL):
#
#   include(joinpath(dirname(dirname(pathof(CutAndDisplaceJulia))), "test", "test_PropagationValleyStress.jl"))
#   cd("path/to/output/folder")
#   (K, E, λ, ν, G) = CutAndDisplaceJulia.ElasticConstantsCheck(ShearModulus(2.0e9), PoissonsRatio(0.25))
#   ρrock = 2900; ρfluid = 2600; Δρ = (ρrock - ρfluid)*9.81
#   testPropagationValleyStress(0, ν, G, Δρ, 5e7, 350; ρrock = ρrock, ValleyWidth = 2000.0, ValleyDepth = -1000.0,
#                               XFactor = 1.0, DepthFactor = 0.5, RadiusFactor = 0.5, VolumeFactor = 1.2, Meshing = :graded)
#
# This file includes test_PropagationTest.jl for DykeBuoyancyGradient, DykeCriticalVolume and the file writers.

include(joinpath(@__DIR__, "test_PropagationTest.jl"))

#########################################################
# Valley (Savage 1984) loaded crack: testPropagationValleyStress and its helpers
#########################################################

"""
    ValleySurfaceXZ(a, b; umax = 20.0, n = 40001)

Valley surface (x, z) from the Savage (1984) conformal map z = w + a*b/(w - i*a) with w real.
`a` is the width parameter [m] and `b` the (negative) depth of the valley floor at x = 0 [m]: the surface is
z = a²b/(u²+a²) in conformal coordinates u, so far from the valley z → 0, the original flat surface.
x runs from about -umax*a to umax*a (x increases monotonically provided 1 + b/a > 0).
"""
function ValleySurfaceXZ(a, b; umax = 20.0, n = 40001)
	u = collect(range(-umax*a, umax*a, length = n))
	w = complex.(u, 0.0)
	z = w .+ a*b ./ (w .- im*a)
	return real.(z), imag.(z)
end

"""
    ValleySurfaceHeight(x, a, b)

Height z of the valley surface at horizontal positions `x` (linear interpolation on ValleySurfaceXZ;
the first/last value is used beyond its end, where the surface is flat to within a few metres).
"""
function ValleySurfaceHeight(x, a, b)
	(xs, zs) = ValleySurfaceXZ(a, b)
	out = zeros(length(x))
	for k in eachindex(x)
		j = searchsortedlast(xs, x[k])
		if j < 1
			out[k] = zs[1]
		elseif j >= length(xs)
			out[k] = zs[end]
		else
			out[k] = zs[j] + (zs[j+1] - zs[j]) * (x[k] - xs[j]) / (xs[j+1] - xs[j])
		end
	end
	return out
end

"""
    ValleyStressTensor(x, z, rg, νLoad, a, b)

Plane-strain stresses (tension positive) in the rock below an incised valley, Savage (1984) gravity solution
with no tectonic stress. x is horizontal (across the valley), z vertical, y the out-of-plane direction (along
the valley). rg = ρrock*g. νLoad is the Poisson's ratio of the loading solution (0.5 gives σxx = σzz far from the valley).
The out-of-plane stress is σyy = νLoad*(σxx + σzz). Returns (σxx, σyy, σzz, σxz) vectors.
"""
function ValleyStressTensor(x, z, rg, νLoad, a, b)
	(σxx, σzz, σxz, ~, ~) = CutAndDisplaceJulia.Savage1984_GravityValleyStress_CoordsIn(0.0, rg, νLoad, a, b, vec(collect(float.(x))), vec(collect(float.(z))))
	σxx = real.(σxx); σzz = real.(σzz); σxz = real.(σxz)
	σyy = νLoad .* (σxx .+ σzz)
	return σxx, σyy, σzz, σxz
end

# Draws one stress panel (colour = most compressive in-plane stress, bars = its direction) on `ax` for the window
# [xlo,xhi]x[zlo,zhi]; positions are plotted in units of a. Returns the heatmap.
function ValleyStressPanel!(fig, ax, rowcol, xlo, xhi, zlo, zhi, a, b, rg, νLoad, xs, zs; nx = 121, nz = 91, bar = 6)
	gx = collect(range(xlo, xhi, length = nx)); gz = collect(range(zlo, zhi, length = nz))
	GX = repeat(gx, 1, nz); GZ = repeat(gz', nx, 1)
	(sxx, syy, szz, sxz) = ValleyStressTensor(vec(GX), vec(GZ), rg, νLoad, a, b)
	sxx = reshape(sxx, size(GX)); szz = reshape(szz, size(GX)); sxz = reshape(sxz, size(GX))
	#Most compressive in-plane principal stress (tension positive => the smallest) and its direction
	half = (sxx .+ szz) ./ 2; rad = sqrt.(((sxx .- szz) ./ 2).^2 .+ sxz.^2)
	σc = half .- rad
	vx = sxz; vz = σc .- sxx
	nv = sqrt.(vx.^2 .+ vz.^2)
	small = nv .< 1e-12 .* (abs.(σc) .+ 1.0)
	vx[small] .= ifelse.(sxx[small] .< szz[small], 1.0, 0.0)
	vz[small] .= ifelse.(sxx[small] .< szz[small], 0.0, 1.0)
	nv = sqrt.(vx.^2 .+ vz.^2); vx = vx ./ nv; vz = vz ./ nv
	above = reshape(vec(GZ) .>= ValleySurfaceHeight(vec(GX), a, b), size(GX))
	σcMPa = σc ./ 1e6; σcMPa[above] .= NaN
	hm = Makie.heatmap!(ax, gx ./ a, gz ./ a, σcMPa)
	Makie.Colorbar(fig[rowcol, 2], hm; label = "most compressive in-plane stress [MPa] (tension +)")
	L = 0.8 * bar * (gx[2] - gx[1]) / a
	segs = Makie.Point2f[]
	for ix in 1:bar:nx, iz in 1:bar:nz
		above[ix, iz] && continue
		cx = gx[ix] / a; cz = gz[iz] / a
		push!(segs, Makie.Point2f(cx - L/2*vx[ix, iz], cz - L/2*vz[ix, iz]))
		push!(segs, Makie.Point2f(cx + L/2*vx[ix, iz], cz + L/2*vz[ix, iz]))
	end
	Makie.linesegments!(ax, segs; color = :white, linewidth = 1.5)
	keep = (xs .>= xlo - 0.1*(xhi - xlo)) .& (xs .<= xhi + 0.1*(xhi - xlo))
	Makie.lines!(ax, xs[keep] ./ a, zs[keep] ./ a; color = :black, linewidth = 3)
	return hm
end

"""
    ValleySetupFigure(filename, P1, P2, P3, a, b, rg, νLoad, xc, zc, R, CritRadius, CritVolume, CrackVolume)

One picture of the setup before the loop, positions in units of the valley width a. Top: x–z section (through the
crack centre) showing the valley surface, the crack and the most compressive in-plane stress (colour, MPa) with its
direction (white bars). The window always contains the valley centre and the crack. Bottom: a zoom on the crack
(orange box in the top panel) with the mesh faces edge-on (red), so the crack is visible whatever its size.
"""
function ValleySetupFigure(filename, P1, P2, P3, a, b, rg, νLoad, xc, zc, R, CritRadius, CritVolume, CrackVolume)
	(xs, zs) = ValleySurfaceXZ(a, b)
	#Overview window: valley axis, both flanks and the crack, with a margin of 0.5a
	xlo = min(-2*a, xc) - 0.5*a; xhi = max(2*a, xc) + 0.5*a
	zlo = min(b, zc) - 0.5*a;    zhi = 0.15*a
	#Zoom window: crack centre ± 2.5 R (square)
	zw = 2.5*R
	fig = try
		Makie.Figure(size = (1000, 1300))
	catch
		Makie.Figure()
	end
	ax1 = Makie.Axis(fig[1, 1]; xlabel = "x / a (across valley)", ylabel = "z / a (up)",
	                 title = "Setup, a = $(a) m, b = $(b) m\ncrack at x = $(round(xc/a, digits=2))a, z = $(round(zc/a, digits=2))a, R = $(round(R, digits=1)) m (critical R = $(round(CritRadius, digits=1)) m), V = $(round(CrackVolume/CritVolume, sigdigits=3))×critical",
	                 aspect = Makie.DataAspect())
	ValleyStressPanel!(fig, ax1, 1, xlo, xhi, zlo, zhi, a, b, rg, νLoad, xs, zs)
	#Mesh faces edge-on: every triangle edge projected on the x–z plane
	tsegs1 = Makie.Point2f[]
	for k in 1:size(P1, 1)
		for (A, B) in ((P1, P2), (P2, P3), (P3, P1))
			push!(tsegs1, Makie.Point2f(A[k, 1]/a, A[k, 3]/a)); push!(tsegs1, Makie.Point2f(B[k, 1]/a, B[k, 3]/a))
		end
	end
	Makie.linesegments!(ax1, tsegs1; color = :red, linewidth = 2)
	Makie.scatter!(ax1, [xc/a], [zc/a]; color = :red, markersize = 10)
	#Zoom box
	bx = [xc - zw, xc + zw, xc + zw, xc - zw, xc - zw] ./ a; bz = [zc - zw, zc - zw, zc + zw, zc + zw, zc - zw] ./ a
	Makie.lines!(ax1, bx, bz; color = :orange, linewidth = 2)
	Makie.xlims!(ax1, xlo/a, xhi/a); Makie.ylims!(ax1, zlo/a, zhi/a)

	ax2 = Makie.Axis(fig[2, 1]; xlabel = "x / a", ylabel = "z / a", title = "Zoom on the crack (mesh faces edge-on)",
	                 aspect = Makie.DataAspect())
	ValleyStressPanel!(fig, ax2, 2, xc - zw, xc + zw, zc - zw, zc + zw, a, b, rg, νLoad, xs, zs; nx = 61, nz = 61, bar = 5)
	Makie.linesegments!(ax2, tsegs1; color = :red, linewidth = 1.5)
	Makie.scatter!(ax2, [xc/a], [zc/a]; color = :red, markersize = 8)
	Makie.xlims!(ax2, (xc - zw)/a, (xc + zw)/a); Makie.ylims!(ax2, (zc - zw)/a, (zc + zw)/a)
	Makie.save(filename, fig)
	return filename
end

"""
    testPropagationValleyStress(HSFlag, ν, G, Δρ, KCrit, NoTris; ρrock, ValleyWidth = 2000.0, ValleyDepth = -1000.0,
                                XFactor = 1.0, DepthFactor = 0.5, RadiusFactor = 0.5, VolumeFactor = 1.2, νLoad = 0.5)

Like testProp, but the crack is loaded by the stresses in the rock below an incised valley (Savage et al. 1984,
plane strain, gravity only, no tectonic stress) and it runs until it reaches the valley surface or is trapped.

  ν, G, Δρ, KCrit  as testProp (ν is the elastic Poisson's ratio used by the BEM)
  ρrock            rock density [kg/m³] (required); fluid density is ρrock - Δρ/g
  ValleyWidth      `a`, valley width parameter [m]. ValleyDepth `b`, depth of the valley floor at x = 0 [m],
                   NEGATIVE (the surface is z = a²b/(u²+a²) in the conformal coordinate u; b = floor depth, a = width)
  XFactor          start centre at x = XFactor*a, measured from the valley axis (x = 0) along the slope direction
                   (both distances are in units of the valley width a)
  DepthFactor      start centre DepthFactor*a below the valley floor level (z = b)
  RadiusFactor     start radius = RadiusFactor * critical radius (the flat, pre-valley lithostatic gradient, vertical crack)
  VolumeFactor     fluid volume = VolumeFactor * critical volume (same flat reference)
  νLoad            Poisson's ratio of the loading solution (0.5: σxx = σzz far field; σyy = νLoad*(σxx+σzz))

The start crack is a penny whose plane contains the most compressive principal stress (and the intermediate one):
its normal is the least compressive principal direction at the crack centre. x is across the valley, y along it,
z up. Writes `SetupFigure` (a png) before the loop. Returns (CritVolFlag, maxX, minX, maxY, minY, maxZ, minZ);
CritVolFlag 1 = reached the valley surface, 0 = trapped, 2 = stopped, 4 = too many loops.
"""
function testPropagationValleyStress(HSFlag,ν,G,Δρ,KCrit,NoTris; ρrock=nothing,
                  ValleyWidth=2000.0, ValleyDepth=-1000.0, XFactor=1.0, DepthFactor=0.5,
                  RadiusFactor=0.5, VolumeFactor=1.2, νLoad=0.5, SetupFigure="ValleyStress-Setup.png",
                  Meshing=:uniform, FrontRefinement=2, GradedCoarsening=1.4, Growth=1.3)
Meshing in (:uniform,:graded) || error("testPropagationValleyStress: Meshing must be :uniform or :graded")
g=9.81
ρrock === nothing && error("testPropagationValleyStress: give the rock density, e.g. ...; ρrock = 2900)")
ρfluid=ρrock-Δρ/g   #fluid density [kg/m³]
rg=ρrock*g
va=ValleyWidth; vb=ValleyDepth    #Savage a and b (not called a, b: the plotting block in the loop uses those names)
vb<0 || error("ValleyDepth (b) must be negative: it is the height of the valley floor below the original flat surface (z = 0)")
(1+vb/va)>0 || error("ValleyDepth must satisfy b > -ValleyWidth (otherwise the valley sides overhang)")

println("creating func vars")
@info HSFlag ν G Δρ KCrit NoTris va vb XFactor DepthFactor RadiusFactor VolumeFactor νLoad

#Critical radius and volume from the flat, pre-valley lithostatic gradient (vertical crack, σh = σv):
#the fluid (ρfluid) against the rock (ρrock) gradient, Δρ
ΔρDip=DykeBuoyancyGradient(Δρ,90.0;StressRatio=1.0,ρrock=ρrock,g=g)
CritRadius=((3*sqrt(pi)*KCrit)/(8*abs(ΔρDip)))^(2/3)
CritVolume=DykeCriticalVolume(ν,G,KCrit,Δρ;Dip=90.0,StressRatio=1.0,ρrock=ρrock)
StartRadius=RadiusFactor*CritRadius
CrackVolume=VolumeFactor*CritVolume
RadiusFactor<1 || printstyled("RadiusFactor >= 1: the start crack is not smaller than the critical radius\n",color=:yellow)
printstyled("""
testPropagationValleyStress: valley a = $(va) m, b = $(vb) m; νLoad = $(νLoad)
          critical radius $(round(CritRadius,sigdigits=4)) m (start radius $(round(StartRadius,sigdigits=4)) m), critical volume $(round(CritVolume,sigdigits=4)) m³
          fluid volume $(round(CrackVolume,sigdigits=4)) m³ (ratio to critical $(round(CrackVolume/CritVolume,sigdigits=3)))
""",color=:cyan)

#Where the crack starts: x along the slope, z below the valley floor level
xc=XFactor*va
zc=vb-DepthFactor*va
if zc+StartRadius>=ValleySurfaceHeight([xc],va,vb)[1]
	error("The start crack would touch or cross the valley surface: increase DepthFactor or reduce RadiusFactor")
end

#Stress at the crack centre -> principal stresses (EigCalc3D). With tension positive, the most compressive
#principal stress is the smallest value; the crack plane contains it, so the crack normal is the largest (least compressive).
(σxxc,σyyc,σzzc,σxzc)=ValleyStressTensor([xc],[zc],rg,νLoad,va,vb)
(S1,S2,S3,S1dir,S2dir,S3dir)=CutAndDisplaceJulia.EigCalc3D(σxxc,σyyc,σzzc,[0.0],σxzc,[0.0])
PrinVals=[S1[1],S2[1],S3[1]]
PrinDirs=[S1dir[1,:],S2dir[1,:],S3dir[1,:]]
iLeast=argmax(PrinVals); iMost=argmin(PrinVals)
CrackNormal=Float64.(PrinDirs[iLeast]); CrackNormal=CrackNormal./sqrt(sum(abs2,CrackNormal))
CrackNormal[1]<0 && (CrackNormal=-CrackNormal)    #sign is arbitrary
MostCompressiveDir=Float64.(PrinDirs[iMost])
@info "Principal stresses at the crack centre (MPa)" S_most_compressive=PrinVals[iMost]/1e6 S_intermediate=PrinVals[6-iMost-iLeast]/1e6 S_least_compressive=PrinVals[iLeast]/1e6 CrackNormal MostCompressiveDir

#Load the penny mesh, rotate its normal onto the crack normal, scale and move it
SurfaceDir=CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"CircleMesh_1a_500Faces.off")
(Points,Triangles)=CutAndDisplaceJulia.OFFReader(SurfaceDir)
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
MeshNormal=[mean(FaceNormalVector[:,1]),mean(FaceNormalVector[:,2]),mean(FaceNormalVector[:,3])]
MeshNormal=MeshNormal./sqrt(sum(abs2,MeshNormal))
cosAB=sum(MeshNormal.*CrackNormal)
if cosAB>1-1e-12
	#already aligned
elseif cosAB<-1+1e-12
	Points[:,2:4]=-Points[:,2:4]
else
	(Points[:,2],Points[:,3],Points[:,4])=CutAndDisplaceJulia.RotateObject3DAllignVectors(MeshNormal,CrackNormal,Points[:,2],Points[:,3],Points[:,4],0.0,0.0,0.0)
end
Points[:,2:4]=Points[:,2:4].*StartRadius
Points[:,2]=Points[:,2].+xc
Points[:,4]=Points[:,4].+zc
if any(Points[:,4].>=ValleySurfaceHeight(Points[:,2],va,vb))
	error("The start crack would touch or cross the valley surface: increase DepthFactor or reduce RadiusFactor")
end
(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )

#One picture of the setup, before the loop
try
	ValleySetupFigure(SetupFigure,P1,P2,P3,va,vb,rg,νLoad,xc,zc,StartRadius,CritRadius,CritVolume,CrackVolume)
	println("Setup figure written to $(joinpath(pwd(),SetupFigure))")
catch e
	@warn "Could not draw the setup figure (the run continues)" exception=(e,catch_backtrace())
end

lps=2000;

#Define a number of tris you want the mesh to have
(target_edge_length,max_target_edge_length)=
CutAndDisplaceJulia.GetDesiredEdgeLength(P1,P2,P3,NoTris)

#For the remeshing in CGAL to work correctly 
(Triangles,Points)=CutAndDisplaceJulia.CreateTrianglesPointsFromP1P2P3(P1,P2,P3)
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)



#Drawing with debugger
scene=[];limits=[];
Px=[];Py=[];Pz=[];
EdgePoints=[]
draw=0
if draw==0
	println("Drawing off")
end
saveallthemeshes=0
if saveallthemeshes==0
	println("Saving only selected meshes (debug is off)")
end
UsingPlots=1
if UsingPlots==1
	println("Plottingsimple on")
end
#If you run a few times this will help differntiate the results
RandNum=rand(1:100)
p=0
StallingFracture=0
MidPointLastLoop=[NaN];P1LastLoop=[NaN];P2LastLoop=[NaN];P3LastLoop=[NaN];FaceNormalVectorLastLoop=[NaN];AreaLastLoop=[NaN]
MidPointLastLoopCheck=[NaN];P1LastLoopCheck=[NaN];P2LastLoopCheck=[NaN];P3LastLoopCheck=[NaN];FaceNormalVectorLastLoopCheck=[NaN];AreaLastLoopCheck=[NaN]

#Looping from here on in
for i=1:lps



	p=1

	xmid=mean(Points[:,2])
	ymid=mean(Points[:,3])
	zmid=mean(Points[:,4])
	xrange=maximum(Points[:,2])-minimum(Points[:,2])

	OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"$i-$p-BeforePolygonRemshing-$RandNum")
	
	@info StartRadius CritRadius target_edge_length OutputDirectory
	
	#=
	#Remesh using Polygon method in CGAL:
	(OutputDirectory)=BuildCGAL.PolygonRemeshingCGAL(OutputDirectory,target_edge_length)
	println(OutputDirectory)
	(Points,Triangles)=CutAndDisplaceJulia.OFFReader(OutputDirectory)

	if draw==1
		p+=1
		if saveallthemeshes==1
			OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"$i-$p-AfterPolygonRemeshingCGAL-$RandNum")
		end
		#@bp
		scene=CutAndDisplaceJuliaPlots.DrawMeshMakie(P1,P2,P3)
		scene=CutAndDisplaceJuliaPlots.SetLookDirection(scene,xmid,ymid,zmid,xrange,"y")


		Makie.save("$i-$p-AfterIsoCGALRemeshing-$RandNum.png", scene)
		p+=1

	end

	#Remesh edges
	(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.CleanEdgeTris(Points,Triangles)
	if draw==1
		scene=CutAndDisplaceJuliaPlots.DrawMeshMakie(P1,P2,P3)
		scene=CutAndDisplaceJuliaPlots.SetLookDirection(scene,xmid,ymid,zmid,xrange,"y")
		Makie.save("$i-$p-AfterCleaningEdgeTris-$RandNum.png", scene)
		if saveallthemeshes==1
			OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"$i-$p-MeshAfterEdgeClean-PreIsocelise-$RandNum")
		end
		p+=1
	end
	(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.IsosceliseEdgeTris(MidPoint,P1,P2,P3,Triangles,Points,FaceNormalVector)
	=#

	if Meshing==:uniform
	#Remesh using Polygon method in CGAL:
	println("Remeshing once")
	OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"BeforePolygonRemshing")
	(OutputDirectory)=BuildCGAL.PolygonRemeshingCGAL(OutputDirectory,target_edge_length)
	(Points,Triangles)=CutAndDisplaceJulia.OFFReader(OutputDirectory)
	(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.CleanEdgeTris(Points,Triangles)
	(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.IsosceliseEdgeTris(MidPoint,P1,P2,P3,Triangles,Points,FaceNormalVector)


	#Do Twice:
	println("Remeshing twice")
	OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"BeforePolygonRemshing")
	(OutputDirectory)=BuildCGAL.PolygonRemeshingCGAL(OutputDirectory,target_edge_length)
	(Points,Triangles)=CutAndDisplaceJulia.OFFReader(OutputDirectory)
	(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.CleanEdgeTris(Points,Triangles)
	(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.IsosceliseEdgeTris(MidPoint,P1,P2,P3,Triangles,Points,FaceNormalVector)


	println("Remeshing 3x")
	OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"BeforePolygonRemshing")
	(OutputDirectory)=BuildCGAL.PolygonRemeshingCGAL(OutputDirectory,target_edge_length)
	(Points,Triangles)=CutAndDisplaceJulia.OFFReader(OutputDirectory)
	(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.CleanEdgeTris(Points,Triangles)
	(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.IsosceliseEdgeTris(MidPoint,P1,P2,P3,Triangles,Points,FaceNormalVector)

	else
	#Graded mesh: tip elements of size GradedCoarsening*h/FrontRefinement on the front, each with its own
	#equilateral rim triangle, growing to GradedCoarsening*h inside. FrontCleanup first does one uniform remesh at h
	#(as above) so the new front from AdvancingFrontCGAL is evenly spaced before it is refined.
	#CleanEdgeTris and IsosceliseEdgeTris are still run, as in the uniform case; on this mesh they
	#change nothing.
	println("Graded remeshing")
	hg=GradedCoarsening*target_edge_length
	(Points,Triangles)=CutAndDisplaceJulia.GradedRemesh(Points,Triangles,hg,hg/FrontRefinement,
	                       BuildCGAL.GradedRemeshingCGAL;
	                       growth=Growth,FrontCleanup=true)
	(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.CleanEdgeTris(Points,Triangles)
	(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.IsosceliseEdgeTris(MidPoint,P1,P2,P3,Triangles,Points,FaceNormalVector)
	println("Graded mesh: $(size(Triangles,1)) triangles")
	end

	#Recompute target_edge_length
	(target_edge_length,max_target_edge_length)=
	CutAndDisplaceJulia.GetDesiredEdgeLength(P1,P2,P3,NoTris)


	if draw==1
		#@bp
		scene=CutAndDisplaceJuliaPlots.DrawMeshMakie(P1,P2,P3)
		scene=CutAndDisplaceJuliaPlots.SetLookDirection(scene,xmid,ymid,zmid,xrange,"y")
		Makie.save("$i-$p-AfterCleaningAndIsolating-$RandNum.png", scene)
		if saveallthemeshes==1
			OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"$i-$p-MeshBeforeSlipCalc-$RandNum")
		end
		p+=1
	end

	#Check if two conditons are filled for n loops then we say the fracture is trapped
	printstyled("OutsideLoop \n",color=:green)
	if i>1
		nlps=15

		if i==2
			#Set the mesh we check to the latest mesh
			MidPointLastLoopCheck=MidPointLastLoop;
			P1LastLoopCheck=P1LastLoop
			P2LastLoopCheck=P2LastLoop
			P3LastLoopCheck=P3LastLoop
			FaceNormalVectorLastLoopCheck=FaceNormalVectorLastLoop
			( AreaLastLoopCheck,~ ) = CutAndDisplaceJulia.AreaOfTriangle3D(  P1LastLoopCheck[:,1],P1LastLoopCheck[:,2],P1LastLoopCheck[:,3],
																	P2LastLoopCheck[:,1],P2LastLoopCheck[:,2],P2LastLoopCheck[:,3],
																	P3LastLoopCheck[:,1],P3LastLoopCheck[:,2],P3LastLoopCheck[:,3] );
			AreaLastLoopCheck=sum(AreaLastLoopCheck)
		end

		InsidePrevious=CutAndDisplaceJulia.CheckIfInsidePreviousBoundary(MidPointLastLoopCheck,P1LastLoopCheck,P2LastLoopCheck,P3LastLoopCheck,
																		FaceNormalVectorLastLoopCheck,
																		P1,P2,P3,MidPoint,max_target_edge_length)
		printstyled("InsideLoop \n",color=:red)
		if InsidePrevious==true
			StallingFracture+=1;
			printstyled("Added one \n",color=:blue)
		else
			StallingFracture=0
			printstyled("Dropped back to 0 \n",color=:blue)

			#Set the mesh we check to the latest mesh
			MidPointLastLoopCheck=MidPointLastLoop;
			P1LastLoopCheck=P1LastLoop
			P2LastLoopCheck=P2LastLoop
			P3LastLoopCheck=P3LastLoop
			FaceNormalVectorLastLoopCheck=FaceNormalVectorLastLoop
			( AreaLastLoopCheck,~ ) = CutAndDisplaceJulia.AreaOfTriangle3D(  P1LastLoopCheck[:,1],P1LastLoopCheck[:,2],P1LastLoopCheck[:,3],
																	P2LastLoopCheck[:,1],P2LastLoopCheck[:,2],P2LastLoopCheck[:,3],
																	P3LastLoopCheck[:,1],P3LastLoopCheck[:,2],P3LastLoopCheck[:,3] );
			AreaLastLoopCheck=sum(AreaLastLoopCheck)

		end

		println("Current number of stalled loops:")
		println(StallingFracture)
		if StallingFracture>nlps
			#Last little check
			( AreaCurrent,~ ) = CutAndDisplaceJulia.AreaOfTriangle3D(  P1[:,1],P1[:,2],P1[:,3],P2[:,1],P2[:,2],P2[:,3],P3[:,1],P3[:,2],P3[:,3] );
			AreaCurrent=sum(AreaCurrent)
			if AreaCurrent>(AreaLastLoopCheck*2)
				#reset
				StallingFracture=0
				#Set the mesh we check to the latest mesh
				MidPointLastLoopCheck=MidPointLastLoop;
				P1LastLoopCheck=P1LastLoop
				P2LastLoopCheck=P2LastLoop
				P3LastLoopCheck=P3LastLoop
				FaceNormalVectorLastLoopCheck=FaceNormalVectorLastLoop
				( AreaLastLoopCheck,~ ) = CutAndDisplaceJulia.AreaOfTriangle3D(  P1LastLoopCheck[:,1],P1LastLoopCheck[:,2],P1LastLoopCheck[:,3],
																		P2LastLoopCheck[:,1],P2LastLoopCheck[:,2],P2LastLoopCheck[:,3],
																		P3LastLoopCheck[:,1],P3LastLoopCheck[:,2],P3LastLoopCheck[:,3] );
				AreaLastLoopCheck=sum(AreaLastLoopCheck)
				continue
			else

			printstyled("Fracture has stalled \n",color=:green)
			
			(Tris,Pnts)=CutAndDisplaceJulia.CreateTrianglesPointsFromP1P2P3(P1,P2,P3)
			maxX=maximum(Pnts[:,2])
			minX=minimum(Pnts[:,2])
			maxY=maximum(Pnts[:,3])
			minY=minimum(Pnts[:,3])
			maxZ=maximum(Pnts[:,4])
			minZ=minimum(Pnts[:,4])
			CritVolFlag=0# - Trapped
			return(CritVolFlag,maxX,minX,maxY,minY,maxZ,minZ)
			
			break

			end
		end
	end


	n=length(P1[:,1]);
	n2=length(Points[:,1]);

	X=MidPoint[:,1];
	Z=MidPoint[:,3];

	#Boundary conditions (tension positive; z is negative below the original flat surface at z = 0)
	#Stresses: the rock, from the Savage (1984) valley solution at each element centre (x across the valley, z up),
	#plane strain: σyy = νLoad*(σxx+σzz) is the out-of-plane stress. No tectonic stress. The valley surface is free of traction.
	(σxx,σyy,σzz,σxz)=ValleyStressTensor(X,Z,rg,νLoad,va,vb)
	σxy = zeros(n)
	σyz = zeros(n)
	#Tractions: the fluid inside the crack, pushing the walls apart: its weight (fluid pressure
	#increasing with depth). SlipCalculator3D adds the uniform fluid pressure that holds CrackVolume
	#(solved for, so no offset or shift is needed here).
	Tn  = -ρfluid*g.*Z   #fluid weight below the original flat surface (z = 0), the datum of the Savage stresses
	Tss = zeros(n)
	Tds = zeros(n)

	Stress=Stresses(σxx,σyy,σzz,σxy,σxz,σyz);
	Traction=Tractions(Tn,Tss,Tds);
	BoundaryConditions=MixedBoundaryConditions(Stress,Traction)
	BoundaryConditions=MixedBoundaryConditionsFluidVolume(BoundaryConditions,CrackVolume)

	#All a fracture
	FractureElements=fill(1,n)

	#Check if fracture has reached the valley surface
	#if HSFlag==1
		if any([P1[:,3];P2[:,3];P3[:,3]].>ValleySurfaceHeight([P1[:,1];P2[:,1];P3[:,1]],va,vb))
			printstyled("Fracture has hit the valley surface \n",color=:green)

			(Tris,Pnts)=CutAndDisplaceJulia.CreateTrianglesPointsFromP1P2P3(P1,P2,P3)
			maxX=maximum(Pnts[:,2])
			minX=minimum(Pnts[:,2])
			maxY=maximum(Pnts[:,3])
			minY=minimum(Pnts[:,3])
			maxZ=maximum(Pnts[:,4])
			minZ=minimum(Pnts[:,4])
			CritVolFlag=1# - Made it to the surface
			return(CritVolFlag,maxX,minX,maxY,minY,maxZ,minZ)

			break
		end
	#end


	@info HSFlag ν G Δρ KCrit CrackVolume NoTris
	Gg=ShearModulus(G); 
	νv=PoissonsRatio(ν);#println("Pr is close to 0.5")
	(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(Gg,νv);

	#try	
	#Calculate slip on faces
	(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FractureElements);

	@info HSFlag ν G Δρ KCrit CrackVolume NoTris

	test=1
	@info test maximum(Dn) minimum(Dn) maximum(Dss) maximum(Dds)


	#Drop interpenetrating elements to 0
	Dn[Dn.<0].=0.0 #Neg Els if they exist dropped to 0


	println("Removeme")
	p+=1
	#CutAndDisplaceJulia.ExportCrackMesh(P1,P2,P3,Dn,Dss,Dds,"$i-$p-MeshFilled-$RandNum")

	if draw==1
		CutAndDisplaceJulia.ExportCrackMesh(P1,P2,P3,Dn,Dss,Dds,"$i-$p-MeshFilled-$RandNum")
		p+=1
		
		scene=CutAndDisplaceJuliaPlots.DrawMeshMakieFilledFaces(Dn,"Dn",P1,P2,P3,cmap2)
		Makie.save("$i-$p-SlipDistribution-$RandNum.png", scene)
		p+=1
	end

	#Compute average face normal
	#AvgFaceNormalVector=[mean(FaceNormalVector[:,1]) mean(FaceNormalVector[:,2]) mean(FaceNormalVector[:,3])]

	#Get tip elements
	(FeP1P2S,FeP1P3S,FeP2P3S)=CutAndDisplaceJulia.GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector)
	#Preclosure and playing with mesh - For running around boundary soon
	P1pre=copy(P1)
	P2pre=copy(P2)
	P3pre=copy(P3)
	MidPointpre=copy(MidPoint)
	FaceNormalVectorpre=copy(FaceNormalVector)

	ClosedEls=round.(Dn,digits=14).==0.0 #eps() <- 16
	nonNan=ClosedEls.==false;
	if any(ClosedEls)
		#Check that the holes touch the edge otherwise dont delete these tris
		P1Closed=copy(P1[ClosedEls,:])
		P2Closed=copy(P2[ClosedEls,:])
		P3Closed=copy(P3[ClosedEls,:])
		(Tris,Pnts)=CutAndDisplaceJulia.CreateTrianglesPointsFromP1P2P3(P1Closed,P2Closed,P3Closed)
		OutputDirectory=CutAndDisplaceJulia.OFFExport(Pnts,Tris,length(Tris[:,1]),length(Pnts[:,1]),"$i-$p-CheckingHoleConnections-$RandNum")
		println(OutputDirectory)
		(OutputDirectory)=BuildCGAL.ConnectedComponentsCGAL(OutputDirectory)
		Flags=CutAndDisplaceJulia.ConnectedComponentsReader(OutputDirectory)
		(SortedTriangles,ConnectedEdge)=CutAndDisplaceJulia.ConnectedConstraints(P1,P2,P3,MidPoint)
	    #number of connected tris (Sorted tris rows not == to 0)
	    NoConnections=sum(SortedTriangles.!=0,dims=2)
		NoConnectedComponents=maximum(Flags)
		ClosedElLocs=findall(ClosedEls)
		for i=1:NoConnectedComponents
			CurrentIndx=findall(Flags.==i)
			if all(NoConnections[ClosedElLocs[CurrentIndx,:]].==3)
				#Set so they are not closed
				ClosedEls[ClosedElLocs[CurrentIndx]].=false
			else 
				#Close this part of the existing edge
				continue
			end
		end
		P1[ClosedEls,:].=NaN
		P2[ClosedEls,:].=NaN
		P3[ClosedEls,:].=NaN
		MidPoint[ClosedEls,:].=NaN
		FaceNormalVector[ClosedEls,:].=NaN
		nonNan=ClosedEls.==false;
	end

	( Area,HalfPerimeter ) = CutAndDisplaceJulia.AreaOfTriangle3D( P1[:,1],P1[:,2],P1[:,3],P2[:,1],P2[:,2],P2[:,3],P3[:,1],P3[:,2],P3[:,3] );
	


	#Check that we have seperated into multiple meshes - if we have fix this
	P1Open=copy(P1[nonNan,:])
	P2Open=copy(P2[nonNan,:])
	P3Open=copy(P3[nonNan,:])
	DnOpen=copy(Dn[nonNan,:])
	(Tris,Pnts)=CutAndDisplaceJulia.CreateTrianglesPointsFromP1P2P3(P1Open,P2Open,P3Open)
	OutputDirectory=CutAndDisplaceJulia.OFFExport(Pnts,Tris,length(Tris[:,1]),length(Pnts[:,1]),"$i-$p-CheckingConnectedComps-$RandNum")
	println(OutputDirectory)
	(OutputDirectory)=BuildCGAL.ConnectedComponentsCGAL(OutputDirectory)
	Flags=CutAndDisplaceJulia.ConnectedComponentsReader(OutputDirectory)
	if any(Flags.>1) #more than one component
		NoConnectedComponents=maximum(Flags)

		ComponentAreas=zeros(NoConnectedComponents)
		ComponentVolumes=zeros(NoConnectedComponents)
		( AreaOpen,~ ) = CutAndDisplaceJulia.AreaOfTriangle3D( P1Open[:,1],P1Open[:,2],P1Open[:,3],P2Open[:,1],P2Open[:,2],P2Open[:,3],P3Open[:,1],P3Open[:,2],P3Open[:,3] );
		for i=1:NoConnectedComponents
			CurrentIndx=findall(Flags.==i)
			ComponentAreas[i]=sum(AreaOpen[CurrentIndx,:])
			ComponentVolumes[i]=sum(AreaOpen[CurrentIndx,:].*DnOpen[CurrentIndx,:])
		end
		#The actual fracture (index in Flags), assumed to be the largest area connected component
		(~,goodflag)=findmax(ComponentAreas) 

		if NoConnectedComponents>1
			printstyled("Volume was $CrackVolume \n",color=:orange)
			for i=1:NoConnectedComponents
				if i!=goodflag
					if isnan(ComponentVolumes[i])
						continue #skip as was clearly a dud tri with 0 area
					end
					CrackVolume=CrackVolume-ComponentVolumes[i]
				end
			end
			printstyled("Parts of fracture have been left behind \n",color=:green)
			printstyled("Volume now $CrackVolume \n",color=:green)
			sleep(1)
		end
		#Indexs of split mesh we dont want:
		BadComponents=findall(Flags.!=goodflag)
		#Set these to closed elements too
		nonNanidx=findall(nonNan)
		ClosedEls[nonNanidx[BadComponents]].=true
		P1[ClosedEls,:].=NaN
		P2[ClosedEls,:].=NaN
		P3[ClosedEls,:].=NaN
		MidPoint[ClosedEls,:].=NaN
		FaceNormalVector[ClosedEls,:].=NaN
		nonNan=ClosedEls.==false;


	end


	if draw==1
		(Tris2,Pnts2)=CutAndDisplaceJulia.CreateTrianglesPointsFromP1P2P3(P1[nonNan,:],P2[nonNan,:],P3[nonNan,:])
		if saveallthemeshes==1
			IgnoreMe=CutAndDisplaceJulia.OFFExport(Pnts2,Tris2,length(Tris2[:,1]),length(Pnts2[:,1]),"$i-$p-ClosedElsRemoved-$RandNum")
		end
		p+=1
	end


	#Get tip elements again now we have closed some elements
	(NewFeP1P2,NewFeP1P3,NewFeP2P3)=CutAndDisplaceJulia.GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector)

	#Just saying the edges we removed we will not comp stress intensity on (or propagate)
	FeP1P2S.FreeFlg[ClosedEls].=false
	FeP1P3S.FreeFlg[ClosedEls].=false
	FeP2P3S.FreeFlg[ClosedEls].=false

	#Comp stress intensity on the old set of elements where the closed els are ignored
	(FeP1P2S,FeP1P3S,FeP2P3S)=CutAndDisplaceJulia.StressIntensity3D(Dn,Dss,Dds,G,ν,FaceNormalVector,FeP1P2S,FeP1P3S,FeP2P3S);

	max_target_edge_length=maximum(HalfPerimeter[nonNan])*(2/3)
	AvgTriangleEdgeLength=mean(HalfPerimeter[nonNan])*(2/3)
	#Graded mesh: the many small tip elements would make the mean edge (and so the largest step a
	#tip can take, in PropagateFracture) much smaller. Use the interior size h so a step is the
	#same length as on the uniform mesh.
	Meshing==:graded && (AvgTriangleEdgeLength=target_edge_length)

	#Switch so new free elements are used but old ones no longer free edges
	FeP1P2S.FreeFlg=NewFeP1P2.FreeFlg;
	FeP1P3S.FreeFlg=NewFeP1P3.FreeFlg;
	FeP2P3S.FreeFlg=NewFeP2P3.FreeFlg;

	#Average STRAIN ENERGY ONLY around border
	avgeverynth=3
	FeP1P2S,FeP1P3S,FeP2P3S=CutAndDisplaceJulia.MovingAverageOfStressIntensity(avgeverynth,P1pre,P2pre,P3pre,FaceNormalVectorpre,MidPointpre,FeP1P2S,FeP1P3S,FeP2P3S)

	#Comp propagation on new els - stress intensity is NaN for new edges
	(p1,p2,p3,StillEdge_P1P2,StillEdge_P1P3,StillEdge_P2P3)=CutAndDisplaceJulia.PropagateFracture( FeP1P2S,FeP1P3S,FeP2P3S,FaceNormalVector,G,ν,KCrit,AvgTriangleEdgeLength )

	#Get vector of strain energy
	StrainEnergyV=[FeP1P2S.StrainEnergy[FeP1P2S.FreeFlg,1]
				  FeP1P3S.StrainEnergy[FeP1P3S.FreeFlg,1]
				  FeP2P3S.StrainEnergy[FeP2P3S.FreeFlg,1]]	

	#Grab the current fracture
	MidPointLastLoop=MidPoint[nonNan,:]
	P1LastLoop=P1[nonNan,:]
	P2LastLoop=P2[nonNan,:]
	P3LastLoop=P3[nonNan,:]
	FaceNormalVectorLastLoop=FaceNormalVector[nonNan,:]

	if UsingPlots==1


		a=1; #plot for y
		b=3; #plot for z
		c=2; #plot for z-
		
		XMid=[FeP1P2S.FeMd[FeP1P2S.FreeFlg,1]
			  FeP1P3S.FeMd[FeP1P3S.FreeFlg,1]
			  FeP2P3S.FeMd[FeP2P3S.FreeFlg,1]]
		YMid=[FeP1P2S.FeMd[FeP1P2S.FreeFlg,2]
			  FeP1P3S.FeMd[FeP1P3S.FreeFlg,2]
			  FeP2P3S.FeMd[FeP2P3S.FreeFlg,2]]
	  	ZMid=[FeP1P2S.FeMd[FeP1P2S.FreeFlg,3]
	  		  FeP1P3S.FeMd[FeP1P3S.FreeFlg,3]
	  	  	  FeP2P3S.FeMd[FeP2P3S.FreeFlg,3]]


		TriDispWriter("$i-$p-MeshWithDisp-ForDrawing-$RandNum",P1LastLoop,P2LastLoop,P3LastLoop,Dn[nonNan,:],Dss[nonNan,:],Dds[nonNan,:])
		CutAndDisplaceJulia.ExportCrackMesh(P1LastLoop,P2LastLoop,P3LastLoop,Dn[nonNan,:],Dss[nonNan,:],Dds[nonNan,:],"$i-$p-MeshFilled-ForDrawing-$RandNum")
		EdgePointsWriter("$i-$p-MeshEdges-ForDrawing-$RandNum",XMid,YMid,ZMid,StrainEnergyV,KCrit)



		# Per-step crack figure (CutAndDisplaceJuliaPlots.PlotCrackState, CairoMakie): x–z, y–z and map
		# views of the opened crack, coloured by opening, tip line coloured by KI/Kc (~0.2 s per step)
		CrackFig=CutAndDisplaceJuliaPlots.PlotCrackState(P1LastLoop,P2LastLoop,P3LastLoop,vec(Dn[nonNan]);
			TipXYZ=[XMid YMid ZMid],TipK=StrainEnergyV,Kc=KCrit,
			meta=(step=i,V=sum(Area[nonNan].*Dn[nonNan])))
		CutAndDisplaceJuliaPlots.SaveFig("$i-$p-CrackViews-$RandNum.png",CrackFig)
		
	end

	#New fracture tip points
	Px=[p1[:,1]; p2[:,1]; p3[:,1]];
	Py=[p1[:,2]; p2[:,2]; p3[:,2]];
	Pz=[p1[:,3]; p2[:,3]; p3[:,3]];

	if isempty(Px)
		printstyled("Fracture has stopped propagating \n",color=:green)

		(Tris,Pnts)=CutAndDisplaceJulia.CreateTrianglesPointsFromP1P2P3(P1,P2,P3)
		maxX=maximum(Pnts[:,2])
		minX=minimum(Pnts[:,2])
		maxY=maximum(Pnts[:,3])
		minY=minimum(Pnts[:,3])
		maxZ=maximum(Pnts[:,4])
		minZ=minimum(Pnts[:,4])
		CritVolFlag=2# - Stopped
		return(CritVolFlag,maxX,minX,maxY,minY,maxZ,minZ)

		break

	end
	
	StillEdge_P1P2[ClosedEls].=false
	StillEdge_P1P3[ClosedEls].=false
	StillEdge_P2P3[ClosedEls].=false

	#previous parts of tip
	if any(StillEdge_P1P2)
		Px=[Px;P1[StillEdge_P1P2,1];P2[StillEdge_P1P2,1]];
		Py=[Py;P1[StillEdge_P1P2,2];P2[StillEdge_P1P2,2]];
		Pz=[Pz;P1[StillEdge_P1P2,3];P2[StillEdge_P1P2,3]];
	end
	if any(StillEdge_P1P3)
		Px=[Px;P1[StillEdge_P1P3,1];P3[StillEdge_P1P3,1]];
		Py=[Py;P1[StillEdge_P1P3,2];P3[StillEdge_P1P3,2]];
		Pz=[Pz;P1[StillEdge_P1P3,3];P3[StillEdge_P1P3,3]];
	end
	if any(StillEdge_P2P3)
		Px=[Px;P2[StillEdge_P2P3,1];P3[StillEdge_P2P3,1]];
		Py=[Py;P2[StillEdge_P2P3,2];P3[StillEdge_P2P3,2]];
		Pz=[Pz;P2[StillEdge_P2P3,3];P3[StillEdge_P2P3,3]];
	end	
	
	EdgePoints=[Px Py Pz]


	OutputDirectory=CutAndDisplaceJulia.xyzExport(EdgePoints[:,1],EdgePoints[:,2],EdgePoints[:,3],"FractureTipPoints")

	if draw==1
		scene=CutAndDisplaceJuliaPlots.DrawMeshMakie(P1[nonNan,:],P2[nonNan,:],P3[nonNan,:])
		scatter!(scene,[Px Py Pz],markersize = 50,limits=scene.limits)#
		scene=CutAndDisplaceJuliaPlots.SetLookDirection(scene,xmid,ymid,zmid,xrange,"y")
		Makie.save("$i-$p-AfterElsRemoved-$RandNum.png", scene)
		p+=1
	end


	#Add to the new points
	Px=[Px; P1[nonNan,1]; P2[nonNan,1]; P3[nonNan,1]];
	Py=[Py; P1[nonNan,2]; P2[nonNan,2]; P3[nonNan,2]];
	Pz=[Pz; P1[nonNan,3]; P2[nonNan,3]; P3[nonNan,3]];

	#remove duplicate rows
	Pnts=[Px Py Pz]
	Pnts=unique(Pnts,dims=1) 

	OutputDirectory=CutAndDisplaceJulia.xyzExport(Pnts[:,1],Pnts[:,2],Pnts[:,3],"NewFracturePoints")
	println(OutputDirectory)
	(OutputDirectory)=BuildCGAL.AdvancingFrontCGAL(OutputDirectory)
	println(OutputDirectory)
	(Points,Triangles)=CutAndDisplaceJulia.OFFReader(OutputDirectory)


	#Clean up advancing front result
	(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
	(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )
	(P1,P2,P3,Points,Triangles)=CutAndDisplaceJulia.RemoveDodgyNewEdges(P1,P2,P3,Points,Triangles,FaceNormalVector,MidPoint,EdgePoints,max_target_edge_length)

	#Export current mesh
	OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"$i-$p-AfterAdvancingFrontClean-$RandNum")
	#Get parameters
	( Area,HalfPerimeter ) = CutAndDisplaceJulia.AreaOfTriangle3D( P1[:,1],P1[:,2],P1[:,3],P2[:,1],P2[:,2],P2[:,3],P3[:,1],P3[:,2],P3[:,3] );
	
	if draw==1
		scene=CutAndDisplaceJuliaPlots.DrawMeshMakie(P1,P2,P3)
		scene=CutAndDisplaceJuliaPlots.SetLookDirection(scene,xmid,ymid,zmid,xrange,"y")
		Makie.save("$i-$p-AfterAdvancingFrontCGALRemeshingAndCleaning-$RandNum.png", scene)
		p+=1
	end

	#Remove extra compoents if they have appeared
	OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"$i-$p-CheckingConnectedCompsP2-$RandNum")
	(OutputDirectory)=BuildCGAL.ConnectedComponentsCGAL(OutputDirectory)
	Flags=CutAndDisplaceJulia.ConnectedComponentsReader(OutputDirectory)
	(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )
	if any(Flags.>1) #more than one component
		NoConnectedComponents=maximum(Flags)
		Areas=zeros(NoConnectedComponents)
		( Area,HalfPerimeter ) = CutAndDisplaceJulia.AreaOfTriangle3D( P1[:,1],P1[:,2],P1[:,3],P2[:,1],P2[:,2],P2[:,3],P3[:,1],P3[:,2],P3[:,3] );
		for i=1:NoConnectedComponents
			CurrentIndx=findall(Flags.==i)
			Areas[i]=sum(Area[CurrentIndx,:])
		end
		#The actual fracture (index in Flags), assumed to be the largest area connected component
		(~,goodflag)=findmax(Areas) 
		#Indexs of split mesh we dont want:
		GoodComponents=findall(Flags.==goodflag)
		#Set these to closed elements too
		P1=copy(P1[GoodComponents,:])
		P2=copy(P2[GoodComponents,:])
		P3=copy(P3[GoodComponents,:])
		(Triangles,Points)=CutAndDisplaceJulia.CreateTrianglesPointsFromP1P2P3(P1,P2,P3)
	end

	(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)

end

printstyled("OverLooped \n",color=:green)

(Tris,Pnts)=CutAndDisplaceJulia.CreateTrianglesPointsFromP1P2P3(P1,P2,P3)
maxX=maximum(Pnts[:,2])
minX=minimum(Pnts[:,2])
maxY=maximum(Pnts[:,3])
minY=minimum(Pnts[:,3])
maxZ=maximum(Pnts[:,4])
minZ=minimum(Pnts[:,4])
CritVolFlag=4# - Stopped
return(CritVolFlag,maxX,minX,maxY,minY,maxZ,minZ)

end
