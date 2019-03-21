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
HSflag=1; #const 

#Traction vector
Tn=zeros(n,1);
Tds=ones(n,1);
Tss=zeros(n,1);

println("Vars created -> to TD func")

(ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn,
 ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
 ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,
 DnUx,DnUy,DnUz,
 DssUx,DssUy,DssUz,
 DdsUx,DdsUy,DdsUz)=
 CutAndDisplaceJulia.TD(x,y,z,P1,P2,P3,DssVec,DdsVec,DnVec,ν,G,DispFlag,StressFlag,HSflag)
 
println("Out of TD func") 
 
#Converting this to stress tensor influences. 
#Uses Hooke's Law to convert strain to stress. Equation 7.131 and 7.132 in David Pollards Book. 
(SxxDss,SyyDss,SzzDss,SxyDss,SxzDss,SyzDss) = CutAndDisplaceJulia.HookesLaw3dStrain2Stress(ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,λ,G);
(SxxDds,SyyDds,SzzDds,SxyDds,SxzDds,SyzDds) = CutAndDisplaceJulia.HookesLaw3dStrain2Stress(ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,λ,G);
(SxxDn,SyyDn,SzzDn,SxyDn,SxzDn,SyzDn) = CutAndDisplaceJulia.HookesLaw3dStrain2Stress(ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn,λ,G);

#Normal traction
DssTn=CutAndDisplaceJulia.CalculateNormalTraction3d( SxxDss,SyyDss,SzzDss,SxyDss,SxzDss,SyzDss,CosAx,CosAy,CosAz )
DdsTn=CutAndDisplaceJulia.CalculateNormalTraction3d( SxxDds,SyyDds,SzzDds,SxyDds,SxzDds,SyzDds,CosAx,CosAy,CosAz )
DnTn =CutAndDisplaceJulia.CalculateNormalTraction3d( SxxDn,SyyDn,SzzDn,SxyDn,SxzDn,SyzDn,CosAx,CosAy,CosAz )

#Cart components of traction vector
(DssT1x,DssT2y,DssT3z ) =CutAndDisplaceJulia.TractionVectorCartesianComponents3d(SxxDss,SyyDss,SzzDss,SxyDss,SxzDss,SyzDss,CosAx,CosAy,CosAz)
(DdsT1x,DdsT2y,DdsT3z ) =CutAndDisplaceJulia.TractionVectorCartesianComponents3d(SxxDds,SyyDds,SzzDds,SxyDds,SxzDds,SyzDds,CosAx,CosAy,CosAz)
(DnT1x,DnT2y,DnT3z ) =CutAndDisplaceJulia.TractionVectorCartesianComponents3d(SxxDn,SyyDn,SzzDn,SxyDn,SxzDn,SyzDn,CosAx,CosAy,CosAz)


#Calculates the directions of the dipslip and ss directions
(StrikeSlipCosine,DipSlipCosine) = CutAndDisplaceJulia.CalculateSSandDSDirs( CosAx,CosAy,CosAz );

#Calculates the directions of the dipslip and ss directions
#StrikeSlipDisplacement_TractionStrikeSlip
( DssTss ) = CutAndDisplaceJulia.CalculateTractionInChosenDirection3d( DssT1x,DssT2y,DssT3z,CosAx,CosAy,CosAz,StrikeSlipCosine );
#DipSlipDisplacement_TractionStrikeSlip
( DdsTss ) = CutAndDisplaceJulia.CalculateTractionInChosenDirection3d( DdsT1x,DdsT2y,DdsT3z,CosAx,CosAy,CosAz,StrikeSlipCosine );
#TensileDisplacement_TractionStrikeSlip
( DnTss  ) = CutAndDisplaceJulia.CalculateTractionInChosenDirection3d( DnT1x,DnT2y,DnT3z,CosAx,CosAy,CosAz,StrikeSlipCosine );

#StrikeSlipDisplacement_TractionDipSlip
( DssTds ) = CutAndDisplaceJulia.CalculateTractionInChosenDirection3d( DssT1x,DssT2y,DssT3z,CosAx,CosAy,CosAz,DipSlipCosine );
#DipSlipDisplacement_TractionDipSlip
( DdsTds ) = CutAndDisplaceJulia.CalculateTractionInChosenDirection3d( DdsT1x,DdsT2y,DdsT3z,CosAx,CosAy,CosAz,DipSlipCosine );
#TensileDisplacement_TractionDipSlip  
( DnTds  ) = CutAndDisplaceJulia.CalculateTractionInChosenDirection3d( DnT1x,DnT2y,DnT3z,CosAx,CosAy,CosAz,DipSlipCosine );

#Now putting influence matricies inside a predefined structure
StressInfMats = StressInf(DnTn,DnTss,DnTds, DssTn,DssTss,DssTds, DdsTn,DdsTss,DdsTds);
DispInfMats   = DispInf(DnUx,DnUy,DnUz,	DssUx,DssUy,DssUz,	DdsUx,DdsUy,DdsUz);

#Concatenate ready for equation 
Atn  = [-StressInfMats.DnTn  -StressInfMats.DssTn -StressInfMats.DdsTn ];     
Atss = [-StressInfMats.DnTss -StressInfMats.DssTss -StressInfMats.DdsTss];     
Atds = [-StressInfMats.DnTds -StressInfMats.DssTds -StressInfMats.DdsTds];     
A= [Atn;Atss;Atds];  

#Prep traction vector
b=[Tn; Tss; Tds];

#Do Equation
D=A\b;