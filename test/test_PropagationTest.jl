function testProp()
#=

	println("PackageCOMPILESTUFF")
	exepath=raw"C:\Users\timmm\AppData\Local\Julia-1.1.0\bin\julia.exe"; 
	new_image=raw"C:\Users\timmm\.julia\dev\PackageCompiler\sysimg\sys.dll"; 
	run(`$exepath --quiet -J $new_image`)

	using AbstractPlotting
	AbstractPlotting.__init__()
	using CutAndDisplaceJuliaPlots
	using Makie

	using DelimitedFiles
	using CutAndDisplaceJulia
	using Statistics
	using BuildCGAL
	using Debugger
	cd(raw"C:\Users\Berlin\Desktop\MeshProp")
	
	using Plots; plot(rand(10))

=#

#Start creating vars for function: 
println("creating func vars")

#Volume
HeightCrack=5; #50
Radius=50
CrackVolume=(π*(Radius^2))*HeightCrack

#Load triangles and points from file (mesh)
#SurfaceDir=CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"VertPenny-300-EqEdges.stl")
SurfaceDir=CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"SquareGrid.stl")
(Points,Triangles)=CutAndDisplaceJulia.STLReader(SurfaceDir)
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
Beta=0; #45 
BetaFromVert=90-Beta;
#Rotate this (YZ)
(Points[:,3],Points[:,4])=CutAndDisplaceJulia.RotateObject2D!(Points[:,3],Points[:,4],0.0,0.0,cosd(BetaFromVert),sind(BetaFromVert))

Radius=Radius/2 #start smaller...
#Get crack to correct radius
Points[:,2:4]=Points[:,2:4].*Radius;
Points[:,4]=Points[:,4].-10000; #2Km deep
(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )
lps=100;

#Define a number of tris you want the mesh to have
NoTris=300;
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
UsingPlots=1
if UsingPlots==1
	println("Plottingsimple on")
end
#If you run a few times this will help differntiate the results
RandNum=rand(1:100)
p=0
#Looping from here on in
for i=1:lps
	p=1

	xmid=mean(Points[:,2])
	ymid=mean(Points[:,3])
	zmid=mean(Points[:,4])
	xrange=maximum(Points[:,2])-minimum(Points[:,2])

	OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"$i-$p-BeforePolygonRemshing-$RandNum")
	
	#Remesh using Polygon method in CGAL:
	(OutputDirectory)=BuildCGAL.PolygonRemeshingCGAL(OutputDirectory,target_edge_length)
	(Points,Triangles)=CutAndDisplaceJulia.OFFReader(OutputDirectory)

	if draw==1
		p+=1
		OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"$i-$p-AfterPolygonRemeshingCGAL-$RandNum")
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
		OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"$i-$p-MeshAfterEdgeClean-PreIsocelise-$RandNum")
		p+=1
	end
	(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.IsosceliseEdgeTrisNew(MidPoint,P1,P2,P3,Triangles,Points,FaceNormalVector)

	#Recompute target_edge_length
	(target_edge_length,max_target_edge_length)=
	CutAndDisplaceJulia.GetDesiredEdgeLength(P1,P2,P3,NoTris)




	if draw==1
		#@bp
		scene=CutAndDisplaceJuliaPlots.DrawMeshMakie(P1,P2,P3)
		scene=CutAndDisplaceJuliaPlots.SetLookDirection(scene,xmid,ymid,zmid,xrange,"y")
		Makie.save("$i-$p-AfterCleaningAndIsolating-$RandNum.png", scene)
		OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"$i-$p-MeshBeforeSlipCalc-$RandNum")
		p+=1
	end

	n=length(P1[:,1]);
	n2=length(Points[:,1]);


	# Which bits we want to compute
	HSFlag=0; #const
	printstyled("Hs off \n",color=:cyan) 

	#Elastic constants
	G=ShearModulus(2.0e9); 
	ν=PoissonsRatio(0.25);#println("Pr is close to 0.5")
	(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(G,ν);

	#Density
	ρrock=2900;
	ρfluid=2600;
	g=9.81;
	Z=MidPoint[:,3];
	Weight=(ρrock*g).*Z;
	#Set BoundaryConditions
	σxx = Weight.*0.8;
	σyy = Weight.*0.8;
	σzz = Weight;
	σxy = zeros(n);    
	σxz = zeros(n);
	σyz = zeros(n);
	Tn=zeros(n)
	#Bouyancy
	println("Check P and T ppr for this eq")
	ZInterface=Z#.+5000
	for j=1:length(ZInterface)
		Tn[j]=ZInterface[j]*(g*(ρrock-ρfluid)); #
	end
	Tss=zeros(n)
	Tds=zeros(n)
	Stress=Stresses(σxx,σyy,σzz,σxy,σxz,σyz);
	Traction=Tractions(Tn,Tss,Tds);
	BoundaryConditions=MixedBoundaryConditions(Stress,Traction)
	BoundaryConditions=MixedBoundaryConditionsFluidVolume(BoundaryConditions,CrackVolume)

	#All a fracture
	FractureElements=fill(1,n)

	#Check if fracture has reached free surface
	if HSFlag==1
		if any([P1[:,3];P2[:,3];P1[:,3]].>0)
			printstyled("Fracture has hit the free surface \n",color=:green)
			break
		end
	end

	#try	
	#Calculate slip on faces
	(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FractureElements);
	
	test=1
	@info test maximum(Dn) minimum(Dn) 

	#Drop interpenetrating elements to 0
	Dn[Dn.<0].=0.0 #Neg Els if they exist dropped to 0

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
		IgnoreMe=CutAndDisplaceJulia.OFFExport(Pnts2,Tris2,length(Tris2[:,1]),length(Pnts2[:,1]),"$i-$p-ClosedElsRemoved-$RandNum")
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
	KCrit=0.05e7; #[5e7 = 50 MPa √m]

	#Switch so new free elements are used but old ones no longer free edges
	FeP1P2S.FreeFlg=NewFeP1P2.FreeFlg;
	FeP1P3S.FreeFlg=NewFeP1P3.FreeFlg;
	FeP2P3S.FreeFlg=NewFeP2P3.FreeFlg;

	#Comp propagation on new els - stress intensity is NaN for new edges
	(p1,p2,p3,StillEdge_P1P2,StillEdge_P1P3,StillEdge_P2P3)=CutAndDisplaceJulia.PropagateFracture( FeP1P2S,FeP1P3S,FeP2P3S,FaceNormalVector,G,ν,KCrit,AvgTriangleEdgeLength )


	if UsingPlots==1
		a=1; #plot for y
		b=3; #plot for z
		XMid=[FeP1P2S.FeMd[FeP1P2S.FreeFlg,a]
			  FeP1P3S.FeMd[FeP1P3S.FreeFlg,a]
			  FeP2P3S.FeMd[FeP2P3S.FreeFlg,a]]
		YMid=[FeP1P2S.FeMd[FeP1P2S.FreeFlg,b]
			  FeP1P3S.FeMd[FeP1P3S.FreeFlg,b]
			  FeP2P3S.FeMd[FeP2P3S.FreeFlg,b]]
		StrainEnergyV=[FeP1P2S.StrainEnergy[FeP1P2S.FreeFlg,1]
					  FeP1P3S.StrainEnergy[FeP1P3S.FreeFlg,1]
					  FeP2P3S.StrainEnergy[FeP2P3S.FreeFlg,1]]	

		fig = plot()
		scatter!([XMid],[YMid],zcolor=StrainEnergyV./KCrit, m=(:blues), lab="")
		display(fig)
		savefig("$i-$p-FaultEdges-$RandNum.png")
		#return XMid, YMid, StrainEnergyV
		
	end

	#New fracture tip points
	Px=[p1[:,1]; p2[:,1]; p3[:,1]];
	Py=[p1[:,2]; p2[:,2]; p3[:,2]];
	Pz=[p1[:,3]; p2[:,3]; p3[:,3]];

	if isempty(Px)
		printstyled("Fracture has stopped propagating \n",color=:green)
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

return P1,P2,P3,Px,Py,Pz

end