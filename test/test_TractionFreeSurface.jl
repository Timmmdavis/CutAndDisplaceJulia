#Test case comparing to Penny shaped crack
function TestBoundaryConditionIsSatisfied(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)
	
	(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);
	 X=MidPoint[:,1];
	 Y=MidPoint[:,2];
	 Z=MidPoint[:,3];



	DispFlag=0;
	StrainFlag=1;

	StrainInfVector=Strains([],[],[],[],[],[]);
	DispInfVector=Disps([],[],[]);
	(StrainAtInfPoints,DispAtInfPoints)= 
	CutAndDisplaceJulia.TD(X,Y,Z,P1,P2,P3,Dss,Dds,Dn,ν,G,DispFlag,StrainFlag,HSFlag,StrainInfVector,DispInfVector)

	#Converting strains to stress tensor influences  
	(Stress) = CutAndDisplaceJulia.HookesLaw3DStrain2Stress(StrainAtInfPoints,λ,G);

	#Extract remote stresses
	if typeof(BoundaryConditions)==MixedBoundaryConditionsFriction
		BoundaryConditions=BoundaryConditions.MixedBoundaryConditions.Stresses
	end
	σxx∞=BoundaryConditions.σxx;
	σyy∞=BoundaryConditions.σyy;
	σzz∞=BoundaryConditions.σzz;
	σxy∞=BoundaryConditions.σxy;
	σxz∞=BoundaryConditions.σxz;
	σyz∞=BoundaryConditions.σyz;
	#Add these to vector
	σxx=Stress.σxx.+σxx∞[1];
	σyy=Stress.σyy.+σyy∞[1];
	σzz=Stress.σzz.+σzz∞[1];
	σxy=Stress.σxy.+σxy∞[1];
	σxz=Stress.σxz.+σxz∞[1];
	σyz=Stress.σyz.+σyz∞[1];

	( Tn,Tds,Tss) = CutAndDisplaceJulia.CalculateNormalAndShearTractions3D( σxx,σyy,σzz,σxy,σxz,σyz,FaceNormalVector);

	@info maximum(abs.(Tn))
	@info maximum(abs.(Tss))
	@info maximum(abs.(Tds))
	if any(abs.([Tn;Tds;Tss]).>1e-8) # Good enough approximation of traction free
		error("The solution is not correctly satisfying constraints that the surface is traction free")
	end	

	#Now testing prop direction -
	#check that new edge points along Md2Ed dir that points along plane we are acting on are twisted right dir
	#Get tip elements
	(FeP1P2S,FeP1P3S,FeP2P3S)=CutAndDisplaceJulia.GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector)
	(FeP1P2S,FeP1P3S,FeP2P3S)=CutAndDisplaceJulia.StressIntensity3D(Dn,Dss,Dds,G,ν,FaceNormalVector,FeP1P2S,FeP1P3S,FeP2P3S);

	(Area,HalfPerimeter ) = CutAndDisplaceJulia.AreaOfTriangle3D( P1[:,1],P1[:,2],P1[:,3],P2[:,1],P2[:,2],P2[:,3],P3[:,1],P3[:,2],P3[:,3] );
	AvgTriangleEdgeLength=mean(HalfPerimeter)*(2/3)

	KCrit=0.0; #[units?]
	(p1,p2,p3,Ang1,Ang2,Ang3)=CutAndDisplaceJulia.PropagateFracture( FeP1P2S,FeP1P3S,FeP2P3S,FaceNormalVector,G,ν,KCrit,AvgTriangleEdgeLength )

	Px=[p1[:,1]; p2[:,1]; p3[:,1]];
	Py=[p1[:,2]; p2[:,2]; p3[:,2]];
	Pz=[p1[:,3]; p2[:,3]; p3[:,3]];

	#Plane normal pointing along x axis. On-in conv - Acting on x
	if FaceNormalVector[1,1]==1
		TwistTest(σxx∞,σxy∞,σxz∞,FeP1P2S,FeP1P3S,FeP2P3S,2,3,Px)
	end
	#Plane normal pointing along y axis. On-in conv - Acting on y
	if FaceNormalVector[1,2]==1
		TwistTest(σyy∞,σxy∞,σyz∞,FeP1P2S,FeP1P3S,FeP2P3S,1,3,Py)
	end
	#Plane normal pointing along z axis. On-in conv - Acting on z
	if FaceNormalVector[1,3]==1
		TwistTest(σzz∞,σxz∞,σyz∞,FeP1P2S,FeP1P3S,FeP2P3S,1,2,Pz)
	end
	#Plane normal pointing along x axis. On-in conv - Acting on x
	if FaceNormalVector[1,1]==-1
		TwistTest(σxx∞,σxy∞,σxz∞,FeP1P2S,FeP1P3S,FeP2P3S,2,3,Px)
	end
	#Plane normal pointing along y axis. On-in conv - Acting on y
	if FaceNormalVector[1,2]==-1
		TwistTest(σyy∞,σxy∞,σyz∞,FeP1P2S,FeP1P3S,FeP2P3S,1,3,Py)
	end
	#Plane normal pointing along z axis. On-in conv - Acting on z
	if FaceNormalVector[1,3]==-1
		TwistTest(σzz∞,σxz∞,σyz∞,FeP1P2S,FeP1P3S,FeP2P3S,1,2,Pz)
	end	

end

function TwistTest(σ11∞,σ12∞,σ13∞,FeP1P2S,FeP1P3S,FeP2P3S,sheartensor2,sheartensor3,Point)
	#function that checks the points at the leading and trailing tip of the crack twist correctly 
	#(wing cracks) when it is subject shear stress/normal stress

	#σ11∞,σ12∞,σ13∞ - the three stresses that act on the plane
	#FeP1P2S - edge points
	#sheartensor2,sheartensor3 - the shear comps acting on the plane
	# i.e. if plane normal faces in x then σxY and σxZ are the shear comps
	#(so values are set to 2 and 3 for Y and Z) (cols ordered in xyz[123])
	#Point - if the plane normal faces in x then this is Px

	if σ11∞[1]==1 
		println("σ11∞ twist test")
		if any(round.(Point.*1e-10).!=0)
			error("Plane not propagating flat")
		end
	end
	if σ12∞[1]==1 
		println("σ12∞ twist test")
		XDir=[FeP1P2S.FeM2Ev[FeP1P2S.FreeFlg,sheartensor2]
			  FeP1P3S.FeM2Ev[FeP1P3S.FreeFlg,sheartensor2]
			  FeP2P3S.FeM2Ev[FeP2P3S.FreeFlg,sheartensor2]]
		#First we check edges pointing along pos y direction
		flag=XDir.==1
		if any(Point[flag].>0) #check no parts prop anticlock wise from tip
			error("Plane not propagating correctly")
		end
		#2nd we check edges pointing along neg y direction
		flag=XDir.==-1
		if any(Point[flag].<0) #check no parts prop anticlock wise from tip
			error("Plane not propagating correctly")
		end		
	end
	if σ13∞[1]==1
		println("σ13∞ twist test")
		YDir=[FeP1P2S.FeM2Ev[FeP1P2S.FreeFlg,sheartensor3]
			  FeP1P3S.FeM2Ev[FeP1P3S.FreeFlg,sheartensor3]
			  FeP2P3S.FeM2Ev[FeP2P3S.FreeFlg,sheartensor3]]
		#First we check edges pointing along pos y direction
		flag=YDir.==1
		if any(Point[flag].>0) #check no parts prop anticlock wise from tip
			error("Plane not propagating correctly")
		end
		#2nd we check edges pointing along neg y direction
		flag=YDir.==-1
		if any(Point[flag].<0) #check no parts prop anticlock wise from tip
			error("Plane not propagating correctly")
		end	
	end
end

function RunBoundaryConditionsForPlane(BndConds1,BndConds2,BndConds3,
	P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,FixedEls)

#Set BoundaryConditions
µ=0.0;
Sf  = 0.0; 
Traction=Tractions(zeros(n),zeros(n),zeros(n))

#Just Stresses as Bnd Cond
TestBoundaryConditionIsSatisfied(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BndConds1,FixedEls)
#Try with friction solver
Friction=FrictionParameters(µ,Sf);
BoundaryConditions=MixedBoundaryConditions(BndConds1,Traction)
BoundaryConditions=MixedBoundaryConditionsFriction(BoundaryConditions,Friction);
TestBoundaryConditionIsSatisfied(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)

#Just Stresses as Bnd Cond
TestBoundaryConditionIsSatisfied(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BndConds2,FixedEls)
#Try with friction solver
Friction=FrictionParameters(µ,Sf);
BoundaryConditions=MixedBoundaryConditions(BndConds2,Traction)
BoundaryConditions=MixedBoundaryConditionsFriction(BoundaryConditions,Friction);
TestBoundaryConditionIsSatisfied(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)

#Just Stresses as Bnd Cond
TestBoundaryConditionIsSatisfied(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BndConds3,FixedEls)
#Try with friction solver
Friction=FrictionParameters(µ,Sf);
BoundaryConditions=MixedBoundaryConditions(BndConds3,Traction)
BoundaryConditions=MixedBoundaryConditionsFriction(BoundaryConditions,Friction);
TestBoundaryConditionIsSatisfied(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)

end



function test_TractionFreeSurfaceFracture()

#Load triangles and points from file (mesh) - Flat xy grid in z plane. 
SurfaceDir=CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"SquareGrid.stl")
(Points,Triangles)=CutAndDisplaceJulia.STLReader(SurfaceDir)
#Very forcefully making sure the surface is flat 
Points[:,4].=0.0; 


#Compute triangle properties
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
n=length(Triangles[:,1]);
n2=length(Points[:,1]);
(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )


#Elastic constants
G=ShearModulus(5000.); 
ν=PoissonsRatio(0.25);
(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(G,ν);
# Which bits we want to compute
HSFlag=0; #const 
FixedEls=zeros(n,1);

σxx = zeros(n);	
σzz = zeros(n);	
σyy = zeros(n);	
σxy = zeros(n);    
σxz = zeros(n);
σyz = zeros(n);
#######Sxz##########
BoundaryConditions1=Stresses(σxx,σyy,σzz,σxy,ones(n),σyz);
#######Syz##########
BoundaryConditions2=Stresses(σxx,σyy,σzz,σxy,σxz,ones(n));
#######Szz##########
BoundaryConditions3=Stresses(σxx,σyy,ones(n),σxy,σxz,σyz);

RunBoundaryConditionsForPlane(BoundaryConditions1,BoundaryConditions2,BoundaryConditions3,
	P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,FixedEls)

##Rotating 180 so normal is other direction
Points[:,2]=-Points[:,2]
#Compute triangle properties
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
n=length(Triangles[:,1]);
n2=length(Points[:,1]);
(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )

RunBoundaryConditionsForPlane(BoundaryConditions1,BoundaryConditions2,BoundaryConditions3,
	P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,FixedEls)


#Now make surface vertical (Switch Y and Z) strikes at 90
PointsTmp=Points[:,4];
Points[:,4]=Points[:,3];
Points[:,3]=PointsTmp;
#Compute triangle properties
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
n=length(Triangles[:,1]);
n2=length(Points[:,1]);
(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )

#######σxy##########
BoundaryConditions1=Stresses(σxx,σyy,σzz,ones(n),σxz,σyz);
#######σyz##########
BoundaryConditions2=Stresses(σxx,σyy,σzz,σxy,σxz,ones(n));
#######σyy##########
BoundaryConditions3=Stresses(σxx,ones(n),σzz,σxy,σxz,σyz);

RunBoundaryConditionsForPlane(BoundaryConditions1,BoundaryConditions2,BoundaryConditions3,
	P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,FixedEls)

##Rotating 180 so normal is other direction
Points[:,2]=-Points[:,2]
#Compute triangle properties
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
n=length(Triangles[:,1]);
n2=length(Points[:,1]);
(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )

RunBoundaryConditionsForPlane(BoundaryConditions1,BoundaryConditions2,BoundaryConditions3,
	P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,FixedEls)


#Now make surface vertical (Switch X and Y) strikes at 0
PointsTmp=Points[:,3];
Points[:,3]=Points[:,2];
Points[:,2]=PointsTmp;
#Compute triangle properties
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
n=length(Triangles[:,1]);
n2=length(Points[:,1]);

(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )

#######σxy##########
BoundaryConditions1=Stresses(σxx,σyy,σzz,ones(n),σxz,σyz);
#######σxz##########
BoundaryConditions2=Stresses(σxx,σyy,σzz,σxy,ones(n),σyz);
#######σxx##########
BoundaryConditions3=Stresses(ones(n),σyy,σzz,σxy,σxz,σyz);

RunBoundaryConditionsForPlane(BoundaryConditions1,BoundaryConditions2,BoundaryConditions3,
	P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,FixedEls)

##Rotating 180 so normal is other direction
Points[:,3]=-Points[:,3]
#Compute triangle properties
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
n=length(Triangles[:,1]);
n2=length(Points[:,1]);
(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )

RunBoundaryConditionsForPlane(BoundaryConditions1,BoundaryConditions2,BoundaryConditions3,
	P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,FixedEls)

end




