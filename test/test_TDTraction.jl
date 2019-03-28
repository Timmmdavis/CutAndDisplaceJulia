#Test case comparing to Penny shaped crack

#Start creating vars for function: 
println("creating func vars")

#Load triangles and points from file (mesh)
ModuleDir=pathof(CutAndDisplaceJulia);
ModuleDir=splitdir(ModuleDir); #remove file name
ModuleDir=ModuleDir[1];
ModuleDir=splitdir(ModuleDir); #out of src
ModuleDir=ModuleDir[1];
if Sys.iswindows()
    SurfaceDir=string(ModuleDir,"\\test\\CircleMesh_1a_500Faces.ts")
else
	SurfaceDir=string(ModuleDir,"/test/CircleMesh_1a_500Faces.ts")
end
(Points,Triangles)=CutAndDisplaceJulia.GoCadAsciiReader(SurfaceDir)

#Compute triangle properties
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
n=length(Triangles[:,1]);
n2=length(Points[:,1]);

(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )

#Elastic constants
G=ShearModulus(1.); 
ν=PoissonsRatio(0.25);
(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(G,ν);


#Repeating array if you want to test with more tris (not the solution vs okada, just speed)
DssVec	=ones(n,1); #const 
DdsVec	=ones(n,1); #const 
DnVec	=ones(n,1); #const 

# Start some vectors (spaced points)
CosAx =  FaceNormalVector[:,1];  
CosAy =  FaceNormalVector[:,2];   
CosAz =  FaceNormalVector[:,3];  
x=zeros(n,1);
y=zeros(n,1);
z=zeros(n,1); 
x[:]=MidPoint[:,1]-(CosAx*1e-12);
y[:]=MidPoint[:,2]-(CosAy*1e-12);
z[:]=MidPoint[:,3]-(CosAz*1e-12);

# What bits we want to compute
DispFlag=1; #const 
StressFlag=1; #const 
HSflag=0; #const 

#Traction vector
Tn=zeros(n,1);
Tds=zeros(n,1);
Tss=ones(n,1);

println("Vars created -> to TD func")

(εxxDn,εyyDn,εzzDn,εxyDn,εxzDn,εyzDn,
 εxxDss,εyyDss,εzzDss,εxyDss,εxzDss,εyzDss,
 εxxDds,εyyDds,εzzDds,εxyDds,εxzDds,εyzDds,
 DnUx,DnUy,DnUz,
 DssUx,DssUy,DssUz,
 DdsUx,DdsUy,DdsUz)=
 CutAndDisplaceJulia.TD(x,y,z,P1,P2,P3,DssVec,DdsVec,DnVec,ν,G,DispFlag,StressFlag,HSflag)
 
println("Out of TD func") 

#Converting strains to stress tensor influences  
(σxxDss,σyyDss,σzzDss,σxyDss,σxzDss,σyzDss) = CutAndDisplaceJulia.HookesLaw3dStrain2Stress(εxxDss,εyyDss,εzzDss,εxyDss,εxzDss,εyzDss,λ,G);
(σxxDds,σyyDds,σzzDds,σxyDds,σxzDds,σyzDds) = CutAndDisplaceJulia.HookesLaw3dStrain2Stress(εxxDds,εyyDds,εzzDds,εxyDds,εxzDds,εyzDds,λ,G);
(σxxDn,σyyDn,σzzDn,σxyDn,σxzDn,σyzDn) 		= CutAndDisplaceJulia.HookesLaw3dStrain2Stress(εxxDn,εyyDn,εzzDn,εxyDn,εxzDn,εyzDn,λ,G);

#Compute normal traction
DssTn=CutAndDisplaceJulia.CalculateNormalTraction3d( σxxDss,σyyDss,σzzDss,σxyDss,σxzDss,σyzDss,CosAx,CosAy,CosAz )
DdsTn=CutAndDisplaceJulia.CalculateNormalTraction3d( σxxDds,σyyDds,σzzDds,σxyDds,σxzDds,σyzDds,CosAx,CosAy,CosAz )
DnTn =CutAndDisplaceJulia.CalculateNormalTraction3d( σxxDn,σyyDn,σzzDn,σxyDn,σxzDn,σyzDn,CosAx,CosAy,CosAz )

#Cart components of traction vector
(DssT1x,DssT2y,DssT3z ) =CutAndDisplaceJulia.TractionVectorCartesianComponents3d(σxxDss,σyyDss,σzzDss,σxyDss,σxzDss,σyzDss,CosAx,CosAy,CosAz)
(DdsT1x,DdsT2y,DdsT3z ) =CutAndDisplaceJulia.TractionVectorCartesianComponents3d(σxxDds,σyyDds,σzzDds,σxyDds,σxzDds,σyzDds,CosAx,CosAy,CosAz)
(DnT1x,DnT2y,DnT3z ) 	=CutAndDisplaceJulia.TractionVectorCartesianComponents3d(σxxDn,σyyDn,σzzDn,σxyDn,σxzDn,σyzDn,CosAx,CosAy,CosAz)

#Calculates the directions of the dipslip and ss directions
(StrikeSlipCosine,DipSlipCosine) = CutAndDisplaceJulia.CalculateSSandDSDirs( CosAx,CosAy,CosAz );

#Compute strike slip traction
( DssTss ) = CutAndDisplaceJulia.CalculateTractionInChosenDirection3d( DssT1x,DssT2y,DssT3z,CosAx,CosAy,CosAz,StrikeSlipCosine );
( DdsTss ) = CutAndDisplaceJulia.CalculateTractionInChosenDirection3d( DdsT1x,DdsT2y,DdsT3z,CosAx,CosAy,CosAz,StrikeSlipCosine );
( DnTss  ) = CutAndDisplaceJulia.CalculateTractionInChosenDirection3d( DnT1x,DnT2y,DnT3z,CosAx,CosAy,CosAz,StrikeSlipCosine );

#Compute dip slip traction
( DssTds ) = CutAndDisplaceJulia.CalculateTractionInChosenDirection3d( DssT1x,DssT2y,DssT3z,CosAx,CosAy,CosAz,DipSlipCosine );
( DdsTds ) = CutAndDisplaceJulia.CalculateTractionInChosenDirection3d( DdsT1x,DdsT2y,DdsT3z,CosAx,CosAy,CosAz,DipSlipCosine );
( DnTds  ) = CutAndDisplaceJulia.CalculateTractionInChosenDirection3d( DnT1x,DnT2y,DnT3z,CosAx,CosAy,CosAz,DipSlipCosine );

#Now putting influence matricies inside a predefined structure
StressInfMats = StressInf(	DnTn,DnTss,DnTds, 	DssTn,DssTss,DssTds, 	DdsTn,DdsTss,DdsTds);
DispInfMats   = DispInf(	DnUx,DnUy,DnUz,		DssUx,DssUy,DssUz,		DdsUx,DdsUy,DdsUz);

#Concatenate influence matrix ready for equation 
Atn  = [-StressInfMats.DnTn  -StressInfMats.DssTn -StressInfMats.DdsTn ];     
Atss = [-StressInfMats.DnTss -StressInfMats.DssTss -StressInfMats.DdsTss];     
Atds = [-StressInfMats.DnTds -StressInfMats.DssTds -StressInfMats.DdsTds];     
A= [Atn;Atss;Atds];  

#Prep traction vector
b=[Tn; Tss; Tds];

#Do linear equation
D=A\b;

#Extract arrays
Dn=D[1:n];
Dss=D[n+1:2*n];
Dds=D[n*2+1:3*n];

#Compute total shearing
TotalShearing = sqrt.((Dss).^2 .+(Dds).^2);

#Finding distance from 0,0 for code midpoints for plotting. 
(θ,ρ) = CutAndDisplaceJulia.cart2pol(x,y);

#Compute analytical solution for shearing of cracks walls
a=1.0; #Radius of penny
(us)=CutAndDisplaceJulia.Eshelby1957_PennyCrackSlipProfile(G,ν,0.0,1.0,ρ,a)#Shearing
(ut)=CutAndDisplaceJulia.Eshelby1957_PennyCrackSlipProfile(G,ν,1.0,0.0,ρ,a)#Opening

#Compute the percent error between analytical and numerical
ResidualPercent=CutAndDisplaceJulia.BAsPercentOfA(us,TotalShearing);

#Test this has not changed 
using Statistics
if maximum(ResidualPercent)>0.012 || Statistics.mean(ResidualPercent)>0.0055
	error("Residual is too high")
end

println("Test Passed")