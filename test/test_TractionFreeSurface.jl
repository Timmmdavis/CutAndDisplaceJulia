#Test case comparing to Penny shaped crack
function TestBoundaryConditionIsSatisfied(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)
	
	(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);
	 X=MidPoint[:,1];
	 Y=MidPoint[:,2];
	 Z=MidPoint[:,3];



	DispFlag=0;
	StressFlag=1;
	#Compute stresses
	(εxxDn,εyyDn,εzzDn,εxyDn,εxzDn,εyzDn,
	εxxDss,εyyDss,εzzDss,εxyDss,εxzDss,εyzDss,
	εxxDds,εyyDds,εzzDds,εxyDds,εxzDds,εyzDds,
	DnUx,DnUy,DnUz,
	DssUx,DssUy,DssUz,
	DdsUx,DdsUy,DdsUz)= 
	CutAndDisplaceJulia.TD(X,Y,Z,P1,P2,P3,Dss,Dds,Dn,ν,G,DispFlag,StressFlag,HSFlag);

	(εxx,εyy,εzz,εxy,εxz,εyz,Ux,Uy,Uz)=
	CutAndDisplaceJulia.TD_sum(εxxDn, εyyDn, εzzDn, εxyDn, εxzDn, εyzDn,
						εxxDss,εyyDss,εzzDss,εxyDss,εxzDss,εyzDss,
						εxxDds,εyyDds,εzzDds,εxyDds,εxzDds,εyzDds,
						DnUx,DnUy,DnUz,
						DssUx,DssUy,DssUz,
						DdsUx,DdsUy,DdsUz);

	#Converting strains to stress tensor influences  
	(σxx,σyy,σzz,σxy,σxz,σyz) = CutAndDisplaceJulia.HookesLaw3DStrain2Stress(εxx,εyy,εzz,εxy,εxz,εyz,λ,G);

	#Extract remote stresses
	if typeof(BoundaryConditions)==MixedBoundaryConditionsFriction
		BoundaryConditions=BoundaryConditions.Stresses
	end
	σxx∞=BoundaryConditions.σxx;
	σyy∞=BoundaryConditions.σyy;
	σzz∞=BoundaryConditions.σzz;
	σxy∞=BoundaryConditions.σxy;
	σxz∞=BoundaryConditions.σxz;
	σyz∞=BoundaryConditions.σyz;
	#Add these to vector
	σxx=σxx.+σxx∞[1];
	σyy=σyy.+σyy∞[1];
	σzz=σzz.+σzz∞[1];
	σxy=σxy.+σxy∞[1];
	σxz=σxz.+σxz∞[1];
	σyz=σyz.+σyz∞[1];

	( Tn,Tds,Tss) = CutAndDisplaceJulia.CalculateNormalAndShearTractions3D( σxx,σyy,σzz,σxy,σxz,σyz,FaceNormalVector);

	@info maximum(abs.(Tn))
	@info maximum(abs.(Tss))
	@info maximum(abs.(Tds))
	if any(abs.([Tn;Tds;Tss]).>1e-8) # Good enough approximation of traction free
		error("The solution is not correctly satisfying constraints that the surface is traction free")
	end	

end


#Load triangles and points from file (mesh)
SurfaceDir=CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"SquareGrid.stl")
(Points,Triangles)=CutAndDisplaceJulia.STLReader(SurfaceDir)

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
BoundaryConditions=MixedBoundaryConditionsFriction(BoundaryConditions,Traction,Friction);
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
BoundaryConditions=MixedBoundaryConditionsFriction(BoundaryConditions,Traction,Friction);
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
BoundaryConditions=MixedBoundaryConditionsFriction(BoundaryConditions,Traction,Friction);
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
BoundaryConditions=MixedBoundaryConditionsFriction(BoundaryConditions,Traction,Friction);
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
BoundaryConditions=MixedBoundaryConditionsFriction(BoundaryConditions,Traction,Friction);
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
BoundaryConditions=MixedBoundaryConditionsFriction(BoundaryConditions,Traction,Friction);
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
BoundaryConditions=MixedBoundaryConditionsFriction(BoundaryConditions,Traction,Friction);
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
BoundaryConditions=MixedBoundaryConditionsFriction(BoundaryConditions,Traction,Friction);
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
BoundaryConditions=MixedBoundaryConditionsFriction(BoundaryConditions,Traction,Friction);
#Calculate slip on faces
#(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);
TestBoundaryConditionIsSatisfied(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)
println("Sxx is satisfied Friction")
