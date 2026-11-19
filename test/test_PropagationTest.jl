# Buoyancy-driven ascent of a fluid-filled crack (dyke) with a fixed fluid volume: testProp.
#
# Usage (see also README section 4):
#
#   using CutAndDisplaceJulia, CutAndDisplaceJuliaPlots, BuildCGAL
#   include(joinpath(dirname(dirname(pathof(CutAndDisplaceJulia))), "test", "test_PropagationTest.jl"))
#   cd("path/to/output/folder")                      # results and figures are written here
#
#   (K, E, λ, ν, G) = CutAndDisplaceJulia.ElasticConstantsCheck(ShearModulus(2.0e9), PoissonsRatio(0.25))
#   ρrock = 2900; ρfluid = 2600; g = 9.81
#   Δρ = (ρrock - ρfluid) * g                        # [Pa/m]
#   Kc = 5e7                                         # fracture toughness [Pa√m] (50 MPa√m)
#   Dip = 60.0                                       # dip of the starting crack [deg] (90 = vertical)
#   StressRatio = 1.0                                # horizontal / vertical rock stress
#   CrackVolume = DykeCriticalVolume(ν, G, Kc, Δρ; Dip = Dip, StressRatio = StressRatio, ρrock = ρrock)
#   (PropFlag, maxX, minX, maxY, minY, maxZ, minZ) =
#       testProp(0, ν, G, Δρ, Kc, CrackVolume, 300; ρrock = ρrock, Dip = Dip, StressRatio = StressRatio)
#
# PropFlag (CritVolFlag): 1 = the crack reached the surface, 0 = it was trapped (the run stops).

#########################################################
using Printf
using Statistics
using CutAndDisplaceJulia, CutAndDisplaceJuliaPlots, BuildCGAL
# Makie is installed as part of CutAndDisplaceJuliaPlots (only used when draw=1 below)
import CutAndDisplaceJuliaPlots: Makie

"""
    DykeBuoyancyGradient(Δρ, Dip; StressRatio = 1.0, ρrock = nothing, g = 9.81)

Pressure gradient [Pa/m] driving the crack up-dip: the fluid against the rock normal stress,
resolved along a crack dipping `Dip` degrees (90 = vertical).
  Isotropic rock stress (StressRatio = 1):  Δρ*sin(Dip)
  σh = StressRatio*σv:  g*[ρrock*(cos²Dip + StressRatio*sin²Dip) - ρfluid]*sin(Dip),  ρfluid = ρrock - Δρ/g
"""
function DykeBuoyancyGradient(Δρ, Dip; StressRatio = 1.0, ρrock = nothing, g = 9.81)
    s = sind(Dip)
    if StressRatio == 1
        return Δρ * s
    end
    ρrock === nothing && error("give ρrock when StressRatio != 1")
    ρfluid = ρrock - Δρ / g
    ρeff_normal = ρrock * (cosd(Dip)^2 + StressRatio * s^2)
    return g * (ρeff_normal - ρfluid) * s
end

"""
    DykeCriticalVolume(ν, G, KCrit, Δρ; Dip = 90.0, StressRatio = 1.0, ρrock = nothing)

Fluid volume [m³] at which a buoyant penny-shaped crack starts to ascend (top tip at KCrit), using
the driving gradient along the crack. A dipping crack (or σh < σv) has less drive along it and so
needs a LARGER volume: e.g. a 45° crack in isotropic stress needs about 1.8x the vertical volume.
"""
function DykeCriticalVolume(ν, G, KCrit, Δρ; Dip = 90.0, StressRatio = 1.0, ρrock = nothing)
    Δρd = DykeBuoyancyGradient(Δρ, Dip; StressRatio = StressRatio, ρrock = ρrock)
    Δρd > 0 || error("No net buoyancy along the crack for these settings (gradient = $Δρd Pa/m)")
    return real(((1 - ν) / (16 * G)) * ((9 * π^4 * KCrit^8) / (complex(Δρd)^5))^(1 / 3))
end

"""
    testProp(HSFlag, ν, G, Δρ, KCrit, CrackVolume, NoTris; ρrock, Dip = 90.0, StressRatio = 1.0)

Propagates a fluid-filled crack holding a fixed volume, from a penny-shaped start, until it reaches
the free surface (z = 0) or is trapped (stops growing).

  HSFlag       0 = full space; 1 = half space (free-surface effect on the stresses)
  ν, G         Poisson's ratio and shear modulus [Pa] (from ElasticConstantsCheck)
  Δρ           (ρrock - ρfluid)*g [Pa/m]
  KCrit        fracture toughness [Pa√m]
  CrackVolume  fluid volume [m³]; DykeCriticalVolume gives the volume needed to ascend at a given dip
  NoTris       approximate number of triangles in the crack mesh
  ρrock        rock density [kg/m³] (required); the fluid density is ρrock - Δρ/g
  Dip          dip of the starting penny [deg], 90 = vertical
  StressRatio  horizontal / vertical rock stress (1 = isotropic)

Loading (z up, z < 0 below the surface; everything opening-positive, as in SlipCalculator3D, where a larger
normal load opens the crack more): the rock stress (σzz = ρrock*g*z, σxx = σyy = StressRatio*σzz, tension
positive) resolved on the crack, plus the weight of the fluid inside it as a traction, Tn = -ρfluid*g*z.
For StressRatio = 1 the net load is Δρ*z: it increases upwards, so the crack opens most at the top and
climbs. The uniform fluid pressure that holds CrackVolume is added by SlipCalculator3D (Brent search).
UnitTests/ut_13_BuoyancyLoading.jl checks this loading and the direction of opening.

Each step writes the opened mesh (*-MeshWithDisp-ForDrawing-*.txt), tip points
(*-MeshEdges-ForDrawing-*.txt) and a figure (*-CrackViews-*.png) to the current folder.
Returns (CritVolFlag, maxX, minX, maxY, minY, maxZ, minZ); CritVolFlag 1 = reached the surface.
"""
function testProp(HSFlag,ν,G,Δρ,KCrit,CrackVolume,NoTris; ρrock=nothing, Dip=90.0, StressRatio=1.0)
g=9.81
ρrock === nothing && error("testProp: give the rock density, e.g. testProp(...; ρrock = 2900)")
ρfluid=ρrock-Δρ/g   #fluid density [kg/m³]


#Start creating vars for function: 
println("creating func vars")
@info HSFlag ν G Δρ KCrit CrackVolume NoTris Dip StressRatio



#Load triangles and points from file (mesh)
#SurfaceDir=CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"VertPenny-300-EqEdges.stl")
SurfaceDir=CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"CircleMesh_1a_500Faces.off")
(Points,Triangles)=CutAndDisplaceJulia.OFFReader(SurfaceDir)
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
#Compute average face normal
#AvgFaceNormalVector=[mean(FaceNormalVector[:,1]) mean(FaceNormalVector[:,2]) mean(FaceNormalVector[:,3])]

#=
#Flatten so normals point up
X=Points[:,2];
Points[:,2]=Points[:,4];
Points[:,4]=X;
=#

#Tilt the penny (rotation in the y-z plane) so it dips Dip degrees (90 = vertical)
(Points[:,3],Points[:,4])=CutAndDisplaceJulia.RotateObject2D!(Points[:,3],Points[:,4],0.0,0.0,cosd(Dip),sind(Dip))

# Radius and volume at which the top tip is critically stressed, using the drive along the crack
ΔρDip=DykeBuoyancyGradient(Δρ,Dip;StressRatio=StressRatio,ρrock=ρrock,g=g)
CritRadius=((3*sqrt(pi)*KCrit)/(8*abs(ΔρDip)))^(2/3)
CritVolume=DykeCriticalVolume(ν,G,KCrit,Δρ;Dip=Dip,StressRatio=StressRatio,ρrock=ρrock)
printstyled("""
testProp: dip $(Dip)°, σh/σv = $(StressRatio)
          fluid volume $(round(CrackVolume,sigdigits=4)) m³, critical volume for this dip $(round(CritVolume,sigdigits=4)) m³ (ratio $(round(CrackVolume/CritVolume,sigdigits=3)))
""",color=:cyan)
CrackVolume<CritVolume && printstyled("The volume is below critical: the crack is likely to be trapped\n",color=:yellow)


StartRadius=CritRadius/2.5

println("Start radius:")
println(StartRadius)

#Get crack to correct radius
Points[:,2:4]=Points[:,2:4].*StartRadius;
Points[:,4]=Points[:,4].-abs(CritRadius*5); #deeper than radius - always below surface
(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )
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
	(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.IsosceliseEdgeTrisNew(MidPoint,P1,P2,P3,Triangles,Points,FaceNormalVector)
	=#

	#Remesh using Polygon method in CGAL:
	println("Remeshing once")
	OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"BeforePolygonRemshing")
	(OutputDirectory)=BuildCGAL.PolygonRemeshingCGAL(OutputDirectory,target_edge_length)
	(Points,Triangles)=CutAndDisplaceJulia.OFFReader(OutputDirectory)
	(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.CleanEdgeTris(Points,Triangles)
	(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.IsosceliseEdgeTrisNew(MidPoint,P1,P2,P3,Triangles,Points,FaceNormalVector)


	#Do Twice:
	println("Remeshing twice")
	OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"BeforePolygonRemshing")
	(OutputDirectory)=BuildCGAL.PolygonRemeshingCGAL(OutputDirectory,target_edge_length)
	(Points,Triangles)=CutAndDisplaceJulia.OFFReader(OutputDirectory)
	(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.CleanEdgeTris(Points,Triangles)
	(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.IsosceliseEdgeTrisNew(MidPoint,P1,P2,P3,Triangles,Points,FaceNormalVector)


	println("Remeshing 3x")
	OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"BeforePolygonRemshing")
	(OutputDirectory)=BuildCGAL.PolygonRemeshingCGAL(OutputDirectory,target_edge_length)
	(Points,Triangles)=CutAndDisplaceJulia.OFFReader(OutputDirectory)
	(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.CleanEdgeTris(Points,Triangles)
	(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.IsosceliseEdgeTrisNew(MidPoint,P1,P2,P3,Triangles,Points,FaceNormalVector)


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

	Z=MidPoint[:,3];

	#Boundary conditions (tension positive; z is negative below the free surface at z = 0)
	#Rock stress: vertical = ρrock*g*z, horizontal = StressRatio * vertical.
	σzz = ρrock*g.*Z
	σxx = StressRatio.*σzz
	σyy = StressRatio.*σzz
	σxy = zeros(n)
	σxz = zeros(n)
	σyz = zeros(n)
	#Fluid weight in the crack, as a traction (opening-positive: fluid pressure pushes the walls apart and
	#increases with depth, -ρfluid*g*z). Together with the rock stress above, resolved on the crack by
	#SlipCalculator3D, the net load is (ρrock*StressRatio-ρfluid)*g*z on a vertical crack: Δρ*z if isotropic.
	Tn  = -ρfluid*g.*Z
	Tss = zeros(n)
	Tds = zeros(n)
	#Stops condition where no pressure will cause crack to fly open: the net normal load (fluid + resolved
	#stress) must not be positive anywhere, or the volume at zero pressure is already too big.
	#(Below the surface with StressRatio*ρrock > ρfluid this is already the case and nothing is shifted.)
	(TnStress,~,~)=CutAndDisplaceJulia.CalculateNormalAndShearTractions3D(σxx,σyy,σzz,σxy,σxz,σyz,FaceNormalVector)
	TnNet=Tn.+TnStress
	@info maximum(TnNet) minimum(TnNet)
	if maximum(TnNet)>0
		Tn=Tn.-maximum(TnNet)
	end

	Stress=Stresses(σxx,σyy,σzz,σxy,σxz,σyz);
	Traction=Tractions(Tn,Tss,Tds);
	BoundaryConditions=MixedBoundaryConditions(Stress,Traction)
	BoundaryConditions=MixedBoundaryConditionsFluidVolume(BoundaryConditions,CrackVolume)

	#All a fracture
	FractureElements=fill(1,n)

	#Check if fracture has reached free surface
	#if HSFlag==1
		if any([P1[:,3];P2[:,3];P3[:,3]].>0)
			printstyled("Fracture has hit the free surface \n",color=:green)

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







function GrabPointNew5(PointsIdxList,P1,P2,P3,j)
#Extract the points on the current bit of the edge
Point=[0. 0. 0.]
Indx=0
for k=1:3
    Idx=PointsIdxList[j,k]
    if Idx==0
        continue
    elseif k==1
        Point=P1[Idx,:]
        Indx=Idx
    elseif k==2
        Point=P2[Idx,:]
        Indx=Idx
    elseif k==3
        Point=P3[Idx,:]
        Indx=Idx
    end
end

return Point,Indx
end


function EdgePointsWriter(filename,X,Y,Z,StrainEnergy,KCrit)
#space delimted text file containing the edges of the mesh and the stresses at these. 

#create a new file for writing
io = open("$filename.txt", "w+");
write(io, "X Y Z StrainEnergy KCrit");print(io,"\r\n");

for i=1:length(X)

	print(io,X[i]);print(io," ");
	print(io,Y[i]);print(io," ");
	print(io,Z[i]);print(io," ");
	print(io,StrainEnergy[i]);print(io," ");
	print(io,KCrit);print(io,"\r\n");

end

print(io,"\r\n");
close(io);

#OutputDir (joinpath uses the right separator on Windows, macOS and Linux)
OutputDirectory=joinpath(pwd(),"$filename.off")

return OutputDirectory
end

function TriDispWriter(filename,P1,P2,P3,Dn,Dss,Dds)
#space delimted text file containing the opening at each tri on the mesh

#create a new file for writing
io = open("$filename.txt", "w+");
write(io, "P1x P1y P1z P2x P2y P2z P3x P3y P3z Dn Dss Dds");print(io,"\r\n");

for i=1:length(P1[:,1])


	print(io,P1[i,1]);print(io," ");
	print(io,P1[i,2]);print(io," ");
	print(io,P1[i,3]);print(io," ");
	print(io,P2[i,1]);print(io," ");
	print(io,P2[i,2]);print(io," ");
	print(io,P2[i,3]);print(io," ");
	print(io,P3[i,1]);print(io," ");
	print(io,P3[i,2]);print(io," ");
	print(io,P3[i,3]);print(io," ");
	print(io,Dn[i]);print(io," ");
	print(io,Dss[i]);print(io," ");
	print(io,Dds[i]);print(io,"\r\n");

end

print(io,"\r\n");
close(io);

#OutputDir (joinpath uses the right separator on Windows, macOS and Linux)
OutputDirectory=joinpath(pwd(),"$filename.off")

return OutputDirectory
end
