#Test case comparing to Penny shaped crack

#Start creating vars for function: 
println("creating func vars")

#Load triangles and points from file (mesh)
SurfaceDir=CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"CircleMesh_1a_500Faces.ts")
(Points,Triangles)=CutAndDisplaceJulia.GoCadAsciiReader(SurfaceDir)

##TestWithSlightDip
#Points[:,4]=Points[:,2].*1e-9;

#Compute triangle properties
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
n=length(Triangles[:,1]);
n2=length(Points[:,1]);

(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )

#Elastic constants
G=ShearModulus(1.); 
ν=PoissonsRatio(0.25);
(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(G,ν);


# Which bits we want to compute
HSFlag=0; #const 


#Traction vector
Tn=zeros(n,1);
Tds=zeros(n,1);
Tss=ones(n,1);

FixedEls=zeros(n,1);

Tn=0.0;
Tds=0.0;
Tss=1.0;

#Set BoundaryConditions
BoundaryConditions=Tractions(Tn,Tss,Tds)

#Friction parameters
µ=0.0;
Sf=0.1;

#Calculate slip on faces
(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3DFriction(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls,µ,Sf);
