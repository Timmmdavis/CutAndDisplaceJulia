function testProp()

#Start creating vars for function: 
println("creating func vars")

#Volume
HeightCrack=2;
Radius=1500
Volume=(π*(Radius^2))*HeightCrack

#Load triangles and points from file (mesh)
SurfaceDir=CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"VertPenny-300-EqEdges.stl")
(Points,Triangles)=CutAndDisplaceJulia.STLReader(SurfaceDir)


#Flatten so normals point up
X=Points[:,2];
Points[:,2]=Points[:,4];
Points[:,4]=X;

#Pennys angle away from Z. 
Beta=-75; 
BetaFromVert=90-Beta;
#Rotate this (YZ)
(Points[:,3],Points[:,4])=CutAndDisplaceJulia.RotateObject2D!(Points[:,3],Points[:,4],0.0,0.0,cosd(BetaFromVert),sind(BetaFromVert))

#Get crack to correct radius
Points[:,2:4]=Points[:,2:4].*Radius;
Points[:,4]=Points[:,4].-2000; #2Km deep
(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )
lps=5;

( Area,HalfPerimeter ) = CutAndDisplaceJulia.AreaOfTriangle3D( P1[:,1],P1[:,2],P1[:,3],P2[:,1],P2[:,2],P2[:,3],P3[:,1],P3[:,2],P3[:,3] );
max_target_edge_length=maximum(HalfPerimeter)*(2/3)
target_edge_length=mean(HalfPerimeter)*(2/3)

#Drawing with debugger
scene=[];limits=[];
Px=[];Py=[];Pz=[];
NewEdgePoints=[]
draw=1

#Looping from here on in
for i=1:lps

	Volume=(π*(Radius^2))*HeightCrack

	#Clean up advancing front result
	(P1,P2,P3,Points,Triangles)=CutAndDisplaceJulia.RemoveDodgyNewEdges(P1,P2,P3,NewEdgePoints,max_target_edge_length)

	#Export current mesh
	OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]))
	println(OutputDirectory)
	#Get parameters
	( Area,HalfPerimeter ) = CutAndDisplaceJulia.AreaOfTriangle3D( P1[:,1],P1[:,2],P1[:,3],P2[:,1],P2[:,2],P2[:,3],P3[:,1],P3[:,2],P3[:,3] );
	
	if draw==1
		scene=CutAndDisplaceJuliaPlots.DrawMeshMakie(P1,P2,P3)
		Makie.save("$i-a-AfterAdvancingFrontCGALRemeshingAndCleaning.png", scene)
	end

	#Remesh using Polygon method in CGAL:
	(OutputDirectory)=BuildCGAL.PolygonRemeshingCGAL(OutputDirectory,target_edge_length)

	#Reload	
	(Points,Triangles)=CutAndDisplaceJulia.OFFReader(OutputDirectory)
	(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )
	#Compute triangle properties
	(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)

	if draw==1
		@bp
		scene=CutAndDisplaceJuliaPlots.DrawMeshMakie(P1,P2,P3)
		Makie.save("$i-b-AfterIsoCGALRemeshing.png", scene)
	end

	#Remesh edges
	(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.CleanAndIsosceliseEdgeTris(MidPoint,P1,P2,P3,Triangles,FaceNormalVector)
	#Recompute target_edge_length
	max_target_edge_length=maximum(HalfPerimeter)*(2/3)
	target_edge_length=mean(HalfPerimeter)*(2/3)

	if draw==1
		@bp
		scene=CutAndDisplaceJuliaPlots.DrawMeshMakie(P1,P2,P3)
		Makie.save("$i-c-AfterCleaningAndIsolating.png", scene)
	end

	n=length(Triangles[:,1]);
	n2=length(Points[:,1]);


	# Which bits we want to compute
	HSFlag=1; #const 

	#Elastic constants
	G=ShearModulus(2.0e9); 
	ν=PoissonsRatio(0.25);
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
	#Bouyancy
	Tn=Z.*(g*(ρrock-ρfluid)); println("Check P and T ppr for this eq")
	Tss=zeros(n)
	Tds=zeros(n)
	Stress=Stresses(σxx,σyy,σzz,σxy,σxz,σyz);
	Traction=Tractions(Tn,Tss,Tds);
	BoundaryConditions=MixedBoundaryConditions(Stress,Traction)
	BoundaryConditions=MixedBoundaryConditionsFluidVolume(BoundaryConditions,Volume)

	#All a fracture
	FractureElements=fill(1,n)

	#try	
	#Calculate slip on faces
	(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FractureElements);

	#Get tip elements
	(FeP1P2S,FeP1P3S,FeP2P3S)=CutAndDisplaceJulia.GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector)

	#Remove closed elements
	Dn[Dn.<0].=0.0 #Neg Els dropped to 0
	Precision=eps(typeof(Dn[1]))*5; #rounding to 5x the system precision
	Precision=1.0/Precision; 
	OpenEls=((round.(Dn.*Precision))./Precision).!=0.0 
	P1=P1[OpenEls,:]
	P2=P2[OpenEls,:]
	P3=P3[OpenEls,:]
	#Just saying the edges we removed we will not comp stress intensity on (or propagate)
	ClosedEls=OpenEls.==false;
	FeP1P2S.FreeFlg[ClosedEls].=false
	FeP1P3S.FreeFlg[ClosedEls].=false
	FeP2P3S.FreeFlg[ClosedEls].=false

	#From tip els too: 
	(FeP1P2S,FeP1P3S,FeP2P3S)=CutAndDisplaceJulia.StressIntensity3D(Dn,Dss,Dds,G,ν,FaceNormalVector,FeP1P2S,FeP1P3S,FeP2P3S);


	( Area,HalfPerimeter ) = CutAndDisplaceJulia.AreaOfTriangle3D( P1[:,1],P1[:,2],P1[:,3],P2[:,1],P2[:,2],P2[:,3],P3[:,1],P3[:,2],P3[:,3] );
	AvgTriangleEdgeLength=mean(HalfPerimeter)*(2/3)
	KCrit=5e7; #[50 MPa √m]
	(p1,p2,p3)=CutAndDisplaceJulia.PropagateFracture( FeP1P2S,FeP1P3S,FeP2P3S,FaceNormalVector,G,ν,KCrit,AvgTriangleEdgeLength )



	Px=[p1[:,1]; p2[:,1]; p3[:,1]];
	Py=[p1[:,2]; p2[:,2]; p3[:,2]];
	Pz=[p1[:,3]; p2[:,3]; p3[:,3]];
	NewEdgePoints=[Px Py Pz]
	@info NewEdgePoints
	OutputDirectory=CutAndDisplaceJulia.xyzExport(NewEdgePoints[:,1],NewEdgePoints[:,2],NewEdgePoints[:,3],"FractureTipPoints")

	if draw==1
		scene=CutAndDisplaceJuliaPlots.DrawMeshMakie(P1,P2,P3)
		scatter!(scene,[Px Py Pz],markersize = 50,limits=scene.limits)#
		Makie.save("$i-d-AfterElsRemoved.png", scene)
	end


	#Add to the new points
	Px=[Px; P1[:,1]; P2[:,1]; P3[:,1]];
	Py=[Py; P1[:,2]; P2[:,2]; P3[:,2]];
	Pz=[Pz; P1[:,3]; P2[:,3]; P3[:,3]];

	#remove duplicate rows
	Pnts=[Px Py Pz]
	Pnts=unique(Pnts,dims=1) 

	OutputDirectory=CutAndDisplaceJulia.xyzExport(Pnts[:,1],Pnts[:,2],Pnts[:,3],"NewFracturePoints")
	println(OutputDirectory)
	(OutputDirectory)=BuildCGAL.AdvancingFrontCGAL(OutputDirectory)
	println(OutputDirectory)
	(Points,Triangles)=CutAndDisplaceJulia.OFFReader(OutputDirectory)
	(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )



end

return P1,P2,P3,Px,Py,Pz

end