#Test case comparing to Penny shaped crack

#Start creating vars for function: 
println("creating func vars")

#Volume
HeightCrack=3;
Radius=1500
Volume=[(π*(Radius^2))*HeightCrack (π*(Radius^2))*HeightCrack*2];

#Load triangles and points from file (mesh)
SurfaceDir=CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"VertPenny-300-EqEdges.stl")
(Points,Triangles)=CutAndDisplaceJulia.STLReader(SurfaceDir)

#Flatten so normals point up
X=Points[:,2];
Points[:,2]=Points[:,4];
Points[:,4]=X;

#Pennys angle away from Z. 
Beta=-30; 
BetaFromVert=90-Beta;
#Rotate this (YZ)
(Points[:,3],Points[:,4])=CutAndDisplaceJulia.RotateObject2D!(Points[:,3],Points[:,4],0.0,0.0,cosd(BetaFromVert),sind(BetaFromVert))

#Get crack to correct radius
Points[:,2:4]=Points[:,2:4].*Radius;
Points[:,4]=Points[:,4].-2000; #2Km deep

#add 2nd crack 
Points2=[Points[:,1:2] Points[:,3].+1000 Points[:,4]] 
Triangles2=copy(Triangles)
#Append to the two meshes
(Points,Triangles,SecondMesh) = CutAndDisplaceJulia.DataAppender3D( Points,Points2,Triangles,Triangles);

#Compute triangle properties
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
n=length(Triangles[:,1]);
n2=length(Points[:,1]);

(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )

# Which bits we want to compute
HSFlag=0; #const 

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
FixedEls=SecondMesh.+1; #single fracture

#Calculate slip on faces
(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);

#Get tip elements
(FeP1P2S,FeP1P3S,FeP2P3S)=CutAndDisplaceJulia.GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector)
(FeP1P2S,FeP1P3S,FeP2P3S)=CutAndDisplaceJulia.StressIntensity3D(Dn,Dss,Dds,G,ν,FaceNormalVector,FeP1P2S,FeP1P3S,FeP2P3S);
