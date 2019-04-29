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

	KCrit=0.0; #[units?]
	(p1,p2,p3)=CutAndDisplaceJulia.PropagateFracture( FeP1P2S,FeP1P3S,FeP2P3S,FaceNormalVector,G,ν,KCrit )

	Px=[p1[:,1]; p2[:,1]; p3[:,1]];
	Py=[p1[:,2]; p2[:,2]; p3[:,2]];
	Pz=[p1[:,3]; p2[:,3]; p3[:,3]];

	println("Do for negative values too")

	#Plane pointing along x axis. On-in conv - Acting on x
	if FaceNormalVector[1,1]==1
		if σxx∞[1]==1 
			if any(round.(Px.*1e-10).!=0)
				error("Plane not propagating flat")
			end
		end
		if σxy∞[1]==1 
			println("σxy∞ twist test")
			YDir=[FeP1P2S.FeM2Ev[FeP1P2S.FreeFlg,2]
				  FeP1P3S.FeM2Ev[FeP1P3S.FreeFlg,2]
				  FeP2P3S.FeM2Ev[FeP2P3S.FreeFlg,2]]
			#First we check edges pointing along pos y direction
			flag=YDir.==1
			#because the top tips of this edge dont progate in the expected direction of a 2D fracture:
			#We compute a ratio between those propagating properly and those not. 
			ratio=sum(Px[flag].<0)/sum(Px[flag].>0)
			if ratio<12 #check no parts prop anticlock wise from tip
				error("Plane not propagating correctly")
			end
			#2nd we check edges pointing along neg y direction
			flag=YDir.==-1
			#because the top tips of this edge dont progate in the expected direction of a 2D fracture:
			#We compute a ratio between those propagating properly and those not. 
			

			
			ratio=sum(Px[flag].>0)/sum(Px[flag].<0)
			println("Good and bad point twist:")
			@info sum(Px[flag].>0) sum(Px[flag].<0)
			#poop
			if ratio<12 #check no parts prop anticlock wise from tip
				error("Plane not propagating correctly")
			end			
		end
		if σxz∞[1]==1
			println("σxz∞ twist test")
			ZDir=[FeP1P2S.FeM2Ev[FeP1P2S.FreeFlg,3]
				  FeP1P3S.FeM2Ev[FeP1P3S.FreeFlg,3]
				  FeP2P3S.FeM2Ev[FeP2P3S.FreeFlg,3]]
			#First we check edges pointing along pos y direction
			flag=ZDir.==1
			#because the top tips of this edge dont progate in the expected direction of a 2D fracture:
			#We compute a ratio between those propagating properly and those not. 
			ratio=sum(Px[flag].<0)/sum(Px[flag].>0)
			if ratio<12 #check no parts prop anticlock wise from tip
				error("Plane not propagating correctly")
			end
			#2nd we check edges pointing along neg y direction
			flag=ZDir.==-1
			#because the top tips of this edge dont progate in the expected direction of a 2D fracture:
			#We compute a ratio between those propagating properly and those not. 
			ratio=sum(Px[flag].>0)/sum(Px[flag].<0)
			if ratio<12 #check no parts prop anticlock wise from tip
				error("Plane not propagating correctly")
			end		
		end
	end
	#Plane pointing along y axis. On-in conv - Acting on y
	if FaceNormalVector[1,2]==1
		if σyy∞[1]==1 
			println("σyy∞ twist test")
			if any(round.(Py.*1e-10).!=0)
				error("Plane not propagating flat")
			end
		end
		if σxy∞[1]==1 

		end
		if σyz∞[1]==1
			
		end
	end
	#Plane pointing along z axis. On-in conv - Acting on z
	if FaceNormalVector[1,3]==1
		if σzz∞[1]==1 
			println("σzz∞ twist test")
			if any(round.(Pz.*1e-10).!=0)
				error("Plane not propagating flat")
			end
		end
		if σxz∞[1]==1 

		end
		if σyz∞[1]==1
			
		end
	end

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

#Set BoundaryConditions
µ=0.0;
Sf  = 0.0; 
σxx = zeros(n);	
σzz = zeros(n);	
σyy = zeros(n);	
σxy = zeros(n);    
σxz = zeros(n);
σyz = zeros(n);
Traction=Tractions(zeros(n),zeros(n),zeros(n))


#######Sxz##########

BoundaryConditions=Stresses(σxx,σyy,σzz,σxy,ones(n),σyz);
#Calculate slip on faces
#(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);
TestBoundaryConditionIsSatisfied(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)
println("Sxz is satisfied")

Friction=FrictionParameters(µ,Sf);
BoundaryConditions=MixedBoundaryConditions(BoundaryConditions,Traction)
BoundaryConditions=MixedBoundaryConditionsFriction(BoundaryConditions,Friction);
#Calculate slip on faces
#(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);
TestBoundaryConditionIsSatisfied(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)
println("Sxz is satisfied Friction")

#######Syz##########

BoundaryConditions=Stresses(σxx,σyy,σzz,σxy,σxz,ones(n));
#Calculate slip on faces
#(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);
TestBoundaryConditionIsSatisfied(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)
println("Syz is satisfied")

Friction=FrictionParameters(µ,Sf);
BoundaryConditions=MixedBoundaryConditions(BoundaryConditions,Traction)
BoundaryConditions=MixedBoundaryConditionsFriction(BoundaryConditions,Friction);
#Calculate slip on faces
#(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);
TestBoundaryConditionIsSatisfied(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)
println("Syz is satisfied Friction")

#######Szz##########
BoundaryConditions=Stresses(σxx,σyy,ones(n),σxy,σxz,σyz);
#Calculate slip on faces
#(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);
TestBoundaryConditionIsSatisfied(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)
println("Szz is satisfied")

Friction=FrictionParameters(µ,Sf);
BoundaryConditions=MixedBoundaryConditions(BoundaryConditions,Traction)
BoundaryConditions=MixedBoundaryConditionsFriction(BoundaryConditions,Friction);
#Calculate slip on faces
#(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);
TestBoundaryConditionIsSatisfied(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)
println("Szz is satisfied Friction")

#Now make surface vertical (Switch Y and Z) strikes at 90
PointsTmp=Points[:,4];
Points[:,4]=Points[:,3];
Points[:,3]=PointsTmp;
#Compute triangle properties
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
n=length(Triangles[:,1]);
n2=length(Points[:,1]);

(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )

 #######Sxy##########

BoundaryConditions=Stresses(σxx,σyy,σzz,ones(n),σxz,σyz);
#Calculate slip on faces
#(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);
TestBoundaryConditionIsSatisfied(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)
println("Sxy is satisfied")

Friction=FrictionParameters(µ,Sf);
BoundaryConditions=MixedBoundaryConditions(BoundaryConditions,Traction)
BoundaryConditions=MixedBoundaryConditionsFriction(BoundaryConditions,Friction);
#Calculate slip on faces
#(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);
TestBoundaryConditionIsSatisfied(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)
println("Sxy is satisfied Friction")

#######σyz##########

BoundaryConditions=Stresses(σxx,σyy,σzz,σxy,σxz,ones(n));
#Calculate slip on faces
#(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);
TestBoundaryConditionIsSatisfied(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)
println("Syz is satisfied")

Friction=FrictionParameters(µ,Sf);
BoundaryConditions=MixedBoundaryConditions(BoundaryConditions,Traction)
BoundaryConditions=MixedBoundaryConditionsFriction(BoundaryConditions,Friction);
#Calculate slip on faces
#(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);
TestBoundaryConditionIsSatisfied(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)
println("Syz is satisfied Friction")

#######Syy##########

BoundaryConditions=Stresses(σxx,ones(n),σzz,σxy,σxz,σyz);
#Calculate slip on faces
#(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);
TestBoundaryConditionIsSatisfied(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)
println("Syy is satisfied")

Friction=FrictionParameters(µ,Sf);
BoundaryConditions=MixedBoundaryConditions(BoundaryConditions,Traction)
BoundaryConditions=MixedBoundaryConditionsFriction(BoundaryConditions,Friction);
#Calculate slip on faces
#(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);
TestBoundaryConditionIsSatisfied(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)
println("Syy is satisfied Friction")

#Now make surface vertical (Switch X and Y) strikes at 0
PointsTmp=Points[:,3];
Points[:,3]=Points[:,2];
Points[:,2]=PointsTmp;
#Compute triangle properties
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
n=length(Triangles[:,1]);
n2=length(Points[:,1]);

(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )


 #######Sxy##########

BoundaryConditions=Stresses(σxx,σyy,σzz,ones(n),σxz,σyz);
#Calculate slip on faces
#(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);
TestBoundaryConditionIsSatisfied(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)
println("Sxy is satisfied")

Friction=FrictionParameters(µ,Sf);
BoundaryConditions=MixedBoundaryConditions(BoundaryConditions,Traction)
BoundaryConditions=MixedBoundaryConditionsFriction(BoundaryConditions,Friction);
#Calculate slip on faces
#(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);
TestBoundaryConditionIsSatisfied(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)
println("Sxy is satisfied Friction")

#######σxz##########

BoundaryConditions=Stresses(σxx,σyy,σzz,σxy,ones(n),σyz);
#Calculate slip on faces
#(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);
TestBoundaryConditionIsSatisfied(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)
println("Sxz is satisfied")

Friction=FrictionParameters(µ,Sf);
BoundaryConditions=MixedBoundaryConditions(BoundaryConditions,Traction)
BoundaryConditions=MixedBoundaryConditionsFriction(BoundaryConditions,Friction);
#Calculate slip on faces
#(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);
TestBoundaryConditionIsSatisfied(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)
println("Sxz is satisfied Friction")



#######Sxx##########

BoundaryConditions=Stresses(ones(n),σyy,σzz,σxy,σxz,σyz);
#Calculate slip on faces
#(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls); 
TestBoundaryConditionIsSatisfied(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)
println("Sxx is satisfied")

Friction=FrictionParameters(µ,Sf);
BoundaryConditions=MixedBoundaryConditions(BoundaryConditions,Traction)
BoundaryConditions=MixedBoundaryConditionsFriction(BoundaryConditions,Friction);
#Calculate slip on faces
#(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);
TestBoundaryConditionIsSatisfied(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)
println("Sxx is satisfied Friction")

end