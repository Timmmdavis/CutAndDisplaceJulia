# Gas-filled crack ascending through the crust: a fixed MASS of ideal gas (compressible), with the
# gas pressure, its expansion and its own weight solved at every step (SlipCalculator3D with
# MixedBoundaryConditionsFluidMass). Rock stress: vertical ρrock*g*z, horizontal StressRatio*vertical.
#
# Usage (needs CutAndDisplaceJuliaPlots, BuildCGAL, Plots):
#   using CutAndDisplaceJulia
#   include(joinpath(dirname(dirname(pathof(CutAndDisplaceJulia))), "test", "test_GasCrackAscent.jl"))
#   cd(mkpath(raw"C:\some\output\folder"))
#   RunGasCrackAscent()                       # air, 500 m, vertical start
#   RunGasCrackAscent(Dip = 45.0)             # 45° start: realigns towards vertical (σh < σv)
#   RunGasCrackAscent(M = 0.016, T = 300.0)   # methane
#
# Note: ideal gas only. CO2 below ~800 m is supercritical and needs a real equation of state.
#
# Outputs in the current folder: GasAscentLog.csv (per step: tip depth, mass, volume, pressure,
# mean density), a CrackViews PNG per step, and GasAscent_Volume_Density.png.

using CutAndDisplaceJulia, CutAndDisplaceJuliaPlots, BuildCGAL, DelimitedFiles

"""
RunGasCrackAscent(; M, T, CrackDepth, Dip, StressRatio, ρrock, G, ν, Kc, StartVolume, NoTris,
                    ArrestMassScaleFactor, MaxIters, patm)

M [kg/mol] gas molar mass (air 0.02897), T [K] temperature, CrackDepth [m] top tip depth at the start,
Dip [deg] of the starting penny (90 = vertical), StressRatio = σh/σv, StartVolume [m³] gas volume at
the start (default: the critical volume for buoyant ascent at that dip). The start radius is chosen so
the tip stress intensity of the pressurised penny is about Kc (no artificial over-pressure).
"""
function RunGasCrackAscent(; M = 0.02897, T = 290.0, CrackDepth = 500.0, Dip = 90.0, StressRatio = 0.75,
                             ρrock = 2700.0, G = 30.0e9, ν = 0.25, Kc = 5e7, StartVolume = nothing,
                             NoTris = 300, ArrestMassScaleFactor = 1.1, MaxIters = 2000, patm = 101325.0)
    g = 9.81
    (Kb, E, λ, νv, Gv) = CutAndDisplaceJulia.ElasticConstantsCheck(ShearModulus(G), PoissonsRatio(ν))
    pref = ρrock * g * CrackDepth + patm                         # absolute pressure at the start depth
    ρgas = pref * M / (8.314462618 * T)                          # gas density there
    # buoyancy gradient along the dipping crack: g*[ρrock*(cos²δ + k sin²δ) - ρgas]*sinδ
    Δρ_dip = g * (ρrock * (cosd(Dip)^2 + StressRatio * sind(Dip)^2) - ρgas) * sind(Dip)
    Δρ_dip > 0 || error("no net buoyancy along the crack for these settings")
    Vcrit = real(((1 - νv) / (16 * Gv)) * ((9 * π^4 * Kc^8) / (complex(Δρ_dip)^5))^(1 / 3))
    V = StartVolume === nothing ? Vcrit : StartVolume
    # penny under uniform pressure: K = 2p√(a/π), V = 16(1-ν²)p a³/(3E)  ->  radius with K = Kc
    a = (3 * E * V / (8 * (1 - νv^2) * sqrt(pi) * Kc))^(2 / 5)
    Mass = ρgas * V
    @info "Gas crack ascent" ρgas Vcrit V a Mass Dip StressRatio
    result = testGasAscent(0, νv, Gv, Kc, Mass, NoTris, a, ArrestMassScaleFactor, CrackDepth;
                           Gas = (M = M, T = T, pref = pref, patm = patm), ρrock = ρrock, g = g,
                           StressRatio = StressRatio, Beta = 90.0 - Dip, MaxIters = MaxIters,
                           LogFile = "GasAscentLog.csv")
    # summary: gas volume and mean density against top-tip depth
    (d, h) = readdlm("GasAscentLog.csv", ',', header = true)
    col(name) = Float64.(d[:, findfirst(==(name), vec(String.(h)))])
    depth = -col("TipZ")
    p1 = Plots.plot(depth, col("Volume"), xflip = true, lw = 2, legend = false,
                    xlabel = "Top tip depth [m]", ylabel = "Gas volume [m³]", title = "Volume")
    p2 = Plots.plot(depth, col("MeanDensity"), xflip = true, lw = 2, legend = false,
                    xlabel = "Top tip depth [m]", ylabel = "Mean gas density [kg/m³]", title = "Density")
    Plots.savefig(Plots.plot(p1, p2, layout = (1, 2), size = (1000, 400)), "GasAscent_Volume_Density.png")
    println("Saved ", joinpath(pwd(), "GasAscent_Volume_Density.png"))
    return result
end

#########################################################
using Printf
using Statistics
using Plots

# Propagation of a fixed MASS of ideal gas (compressible) up through the crust.
#   Gas         (M = molar mass [kg/mol], T = temperature [K], pref = reference absolute pressure [Pa],
#                patm = atmospheric pressure [Pa]); the gas pressure and its own weight are solved for.
#   ρrock, g    rock density and gravity; rock stress σzz = ρrock*g*z, σxx = σyy = StressRatio*σzz
#   Beta        tilt of the starting penny away from vertical [deg] (dip = 90 - Beta)
#   CrackMass   gas mass [kg]; multiplied by ArrestMassScaleFactor each time the crack arrests
function testGasAscent(HSFlag,ν,G,KCrit,CrackMass,NoTris,StartRadius,ArrestMassScaleFactor,CrackDepth;
                       Gas, ρrock=2700.0, g=9.81, StressRatio=0.75, Beta=0.0, MaxIters=2000, LogFile="GasAscentLog.csv")
open(LogFile,"w") do io
	println(io,"Iter,n,TipZ,MinZ,Mass,Volume,P,MeanDensity,iterations,converged,nclosed")
end
LastRho=NaN



#Start creating vars for function: 
println("creating func vars")
@info HSFlag ν G KCrit CrackMass NoTris StartRadius ArrestMassScaleFactor CrackDepth Gas



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

#Pennys angle away from Z. 
#Beta (keyword): tilt of the penny away from vertical [deg]; dip = 90-Beta
BetaFromVert=90-Beta;
#Rotate this (YZ)
(Points[:,3],Points[:,4])=CutAndDisplaceJulia.RotateObject2D!(Points[:,3],Points[:,4],0.0,0.0,cosd(BetaFromVert),sind(BetaFromVert))

# Top tip critically stressed

println("Start radius:")
println(StartRadius)

#Get crack to correct radius
Points[:,2:4]=Points[:,2:4].*StartRadius;
#Place the top tip of the crack CrackDepth below the free surface (z=0),
#so when the top tip reaches z=0 the crack has ascended CrackDepth
Points[:,4]=Points[:,4].-(maximum(Points[:,4])+CrackDepth);
(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )
lps=MaxIters;

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
#Number of times the fracture has arrested and had its volume increased
NoArrests=0
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
	
	@info StartRadius target_edge_length OutputDirectory
	
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

			#Arrested: increase the crack volume and keep going from the current mesh
			NoArrests+=1
			CrackMass=CrackMass*ArrestMassScaleFactor
			printstyled("Arrest #$NoArrests: gas mass increased by x$ArrestMassScaleFactor to $CrackMass kg \n",color=:light_magenta)
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

			end
		end
	end


	n=length(P1[:,1]);
	n2=length(Points[:,1]);

	Z=MidPoint[:,3];

	#Set BoundaryConditions
	σxx = zeros(n);  
	σyy = zeros(n);  
	σzz = zeros(n);  
	σxy = zeros(n);    
	σxz = zeros(n);
	σyz = zeros(n);
	#Bouyancy
	Tn=zeros(n)
	Tss=zeros(n)
	Tds=zeros(n)

	#Rock stress: vertical = ρrock*g*z (z negative below the surface, tension positive),
	#horizontal = StressRatio * vertical. The gas pressure and its weight are solved for, so no
	#fluid traction is applied here.
	σzz=ρrock*g.*Z
	σxx=StressRatio.*σzz
	σyy=StressRatio.*σzz

	Stress=Stresses(σxx,σyy,σzz,σxy,σxz,σyz);
	Traction=Tractions(Tn,Tss,Tds);
	BoundaryConditions=MixedBoundaryConditions(Stress,Traction)
	BoundaryConditions=IdealGasFluidMass(BoundaryConditions,CrackMass,Gas.M,Gas.T,Gas.pref;patm=Gas.patm)

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


	@info HSFlag ν G KCrit CrackMass NoTris
	Gg=ShearModulus(G); 
	νv=PoissonsRatio(ν);#println("Pr is close to 0.5")
	(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(Gg,νv);

	#try	
	#Calculate slip on faces
	#Slip, opening and gas state on the current mesh
	(Dn, Dss, Dds, Fluid)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FractureElements);
	LastRho=Fluid.ρ̄[1]
	( AreaNow,_ ) = CutAndDisplaceJulia.AreaOfTriangle3D( P1[:,1],P1[:,2],P1[:,3],P2[:,1],P2[:,2],P2[:,3],P3[:,1],P3[:,2],P3[:,3] );
	open(LogFile,"a") do io
		ZA=[P1[:,3];P2[:,3];P3[:,3]]
		println(io,"$i,$n,$(maximum(ZA)),$(minimum(ZA)),$CrackMass,$(sum(vec(AreaNow).*Dn)),$(Fluid.P[1]),$(Fluid.ρ̄[1]),$(Fluid.iterations),$(Fluid.converged),$(sum(Dn.<=0))")
	end
	printstyled("Gas: P = $(round(Fluid.P[1]/1e6,sigdigits=4)) MPa, mean density = $(round(Fluid.ρ̄[1],sigdigits=4)) kg/m³ \n",color=:cyan)

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
			printstyled("Mass was $CrackMass kg \n",color=:orange)
			for i=1:NoConnectedComponents
				if i!=goodflag
					if isnan(ComponentVolumes[i])
						continue #skip as was clearly a dud tri with 0 area
					end
					CrackMass=CrackMass-ComponentVolumes[i]*LastRho   #gas left behind (mean density)
				end
			end
			printstyled("Parts of fracture have been left behind \n",color=:green)
			printstyled("Mass now $CrackMass kg \n",color=:green)
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



		# Per-step crack figure (CutAndDisplaceJuliaPlots.PlotCrackState, CairoMakie): ~0.2 s vs ~2.8 s for the old Plots figure
		CrackFig=CutAndDisplaceJuliaPlots.PlotCrackState(P1LastLoop,P2LastLoop,P3LastLoop,vec(Dn[nonNan]);
			TipXYZ=[XMid YMid ZMid],TipK=StrainEnergyV,Kc=KCrit,
			meta=(step=i,V=sum(Area[nonNan].*Dn[nonNan])))
		CutAndDisplaceJuliaPlots.SaveFig("$i-$p-CrackViews-$RandNum.png",CrackFig)

		if Sys.islinux()
			#attempt to stop linux error
			GR.inline("png")
		end
		
	end

	#Report current crack state now we have the displacements, just before remeshing
	ZAll=[P1[nonNan,3];P2[nonNan,3];P3[nonNan,3]]
	HalfCrackHeight=(maximum(ZAll)-minimum(ZAll))/2
	OpenedVolume=sum(Area[nonNan].*Dn[nonNan])
	printstyled("Iter $i | half crack height = $(round(HalfCrackHeight,digits=2)) m | top tip depth = $(round(-maximum(ZAll),digits=2)) m | gas mass = $(round(CrackMass,sigdigits=5)) kg | volume (sum A*Dn) = $(round(OpenedVolume,sigdigits=5)) m³ \n",color=:light_magenta)

	#New fracture tip points
	Px=[p1[:,1]; p2[:,1]; p3[:,1]];
	Py=[p1[:,2]; p2[:,2]; p3[:,2]];
	Pz=[p1[:,3]; p2[:,3]; p3[:,3]];

	if isempty(Px)
		printstyled("Fracture has stopped propagating \n",color=:green)

		#Arrested: increase the crack volume and keep going from the current mesh
		NoArrests+=1
		CrackMass=CrackMass*ArrestMassScaleFactor
		printstyled("Arrest #$NoArrests: gas mass increased by x$ArrestMassScaleFactor to $CrackMass kg \n",color=:light_magenta)
		StallingFracture=0
		continue

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

	Printf.@printf(io,"%e",X[i]);print(io," ");
	Printf.@printf(io,"%e",Y[i]);print(io," ");
	Printf.@printf(io,"%e",Z[i]);print(io," ");
	Printf.@printf(io,"%e",StrainEnergy[i]);print(io," ");
	Printf.@printf(io,"%e",KCrit);print(io,"\r\n");

end

print(io,"\r\n");
close(io);

if Sys.islinux()
	#OutputDir
	OutputDirectory=string(pwd(),"/$filename.off")
else
	#OutputDir
	OutputDirectory=string(pwd(),"\\$filename.off")
end

return OutputDirectory
end

function TriDispWriter(filename,P1,P2,P3,Dn,Dss,Dds)
#space delimted text file containing the opening at each tri on the mesh

#create a new file for writing
io = open("$filename.txt", "w+");
write(io, "P1x P1y P1z P2x P2y P2z P3x P3y P3z Dn Dss Dds");print(io,"\r\n");

for i=1:length(P1[:,1])


	Printf.@printf(io,"%e",P1[i,1]);print(io," ");
	Printf.@printf(io,"%e",P1[i,2]);print(io," ");
	Printf.@printf(io,"%e",P1[i,3]);print(io," ");
	Printf.@printf(io,"%e",P2[i,1]);print(io," ");
	Printf.@printf(io,"%e",P2[i,2]);print(io," ");
	Printf.@printf(io,"%e",P2[i,3]);print(io," ");
	Printf.@printf(io,"%e",P3[i,1]);print(io," ");
	Printf.@printf(io,"%e",P3[i,2]);print(io," ");
	Printf.@printf(io,"%e",P3[i,3]);print(io," ");
	Printf.@printf(io,"%e",Dn[i]);print(io," ");
	Printf.@printf(io,"%e",Dss[i]);print(io," ");
	Printf.@printf(io,"%e",Dds[i]);print(io,"\r\n");

end

print(io,"\r\n");
close(io);

if Sys.islinux()
	#OutputDir
	OutputDirectory=string(pwd(),"/$filename.off")
else
	#OutputDir
	OutputDirectory=string(pwd(),"\\$filename.off")
end

return OutputDirectory
end