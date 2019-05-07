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
FractureElements=zeros(size(Triangles,1))

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
n=5;
#Looping from here on in
for i=1:n

	#Export current mesh
	OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]))
	println(OutputDirectory)
	#Remesh using Polygon method in CGAL:
	( Area,HalfPerimeter ) = CutAndDisplaceJulia.AreaOfTriangle3D( P1[:,1],P1[:,2],P1[:,3],P2[:,1],P2[:,2],P2[:,3],P3[:,1],P3[:,2],P3[:,3] );
	target_edge_length=mean(HalfPerimeter)*(2/3)
	(Outputdirectory)=BuildCGAL.PolygonRemeshingCGAL(OutputDirectory,target_edge_length)

	#Reload	
	(Points,Triangles)=CutAndDisplaceJulia.OFFReader(Outputdirectory)
	(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )

	@bp

	#Compute triangle properties
	(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
	n=length(Triangles[:,1]);
	n2=length(Points[:,1]);
	(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )
	(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.CleanAndIsosceliseEdgeTris(MidPoint,P1,P2,P3,Triangles,FaceNormalVector)

	#kill 
	poop

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
	σxx = Weight;
	σyy = Weight;
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
	FractureElements=FractureElements.+1; #single fracture


	#Calculate slip on faces
	(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FractureElements);

	#Get tip elements
	(FeP1P2S,FeP1P3S,FeP2P3S)=CutAndDisplaceJulia.GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector)
	(FeP1P2S,FeP1P3S,FeP2P3S)=CutAndDisplaceJulia.StressIntensity3D(Dn,Dss,Dds,G,ν,FaceNormalVector,FeP1P2S,FeP1P3S,FeP2P3S);

	KCrit=1e3; #[units?]
	(p1,p2,p3)=CutAndDisplaceJulia.PropagateFracture( FeP1P2S,FeP1P3S,FeP2P3S,FaceNormalVector,G,ν,KCrit )

	Px=[p1[:,1]; p2[:,1]; p3[:,1]];
	Py=[p1[:,2]; p2[:,2]; p3[:,2]];
	Pz=[p1[:,3]; p2[:,3]; p3[:,3]];
	#scatter(P1[:,1],P1[:,2])

end

end