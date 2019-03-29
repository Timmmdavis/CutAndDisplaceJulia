#Test case comparing to Penny shaped crack

#Start creating vars for function: 
println("creating func vars")

#Load triangles and points from file (mesh)
SurfaceDir=CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"ValleySurface_ProperEq_1000Faces.ts")
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


# Which bits we want to compute
HSFlag=0; #const 

#Locked Els
FixedEls= zeros(size(Triangles[:,1])); 

#TEST locking edges MATLAB CODE FOR MOMENT
#TR = triangulation(Triangles,Points(:,2:4));
#[TotalConnectionDistance,SortedTriangles,Noconnections] = ConnectedTrianglesFinder(TR,MidPoint);
#Flag= Noconnections(:,2)<3;
#Fdisp(Flag)=1;


#Set BoundaryConditions
P=1; 		#Density of elastic material 
g=1;		#Acceleration due to gravity
Tectsx=0; 	#Tectonic stress (xx)
σxx = ((ν/(1-ν))*P*g*.+MidPoint[:,3]).+Tectsx;	#Equation 1b, Martel and Muller
σzz = P*g*.+MidPoint[:,3];			#Equation 1a, Martel and Muller
σyy = ν*(σxx.+σzz);					#Plane strain conditions
σxy = 0;          					#Positive is right lateral movement
σxz = 0;
σyz = 0;
BoundaryConditions=Stresses(σxx,σyy,σzz,σxy,σxz,σyz);

#Calculate slip on faces
(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);

(FixedEls,Triangles,FaceNormalVector,MidPoint,P1,P2,P3)=
CutAndDisplaceJulia.RemoveFixedElements3D(FixedEls,Triangles,FaceNormalVector,MidPoint,P1,P2,P3)

#A and B describing slope, See paper for diagram
a=2;
b=-1;
#Setting up grid points U and V, these are mapped to a different location
#later
# Start some vectors (spaced points)
xx = range(0,stop=4,length=50); #linspace deprecated
yy = range(-0.3265,stop=-4,length=46); #linspace deprecated
(u,v)=CutAndDisplaceJulia.meshgrid(xx,yy);
#Analytical solution
(σxxAn,σyyAn,σxyAn,xobs,yobs) = CutAndDisplaceJulia.Savage1984_GravityValleyStress(Tectsx,P*g,ν,a,b,u,v);


X=xobs[:].+50;
Y=zeros(size(X));
Z=yobs[:];
DispFlag=1;
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


#Adding vertical Stress profile: 
σxxGrav = (ν/(1-ν))*P*g.*Z;	#Equation 1b, Martel and Muller
σzzGrav = P*g.*Z;				#Equation 1a, Martel and Muller
σxzGrav =zeros(size(σxxGrav));
σxx=σxxGrav.+σxx;
σzz=σzzGrav.+σzz;
σxz=σxzGrav.+σxz;


#Printing results (residual) - NOTE THE DIFFERENCT COORDINATE SYSTEMS 2D vs 3D
σxxRes=σxx.-σxxAn[:];
σxxRes=σxxRes[isnan.(σxxRes).==0];#removing nans to get proper mean
σzzRes=σzz.-σyyAn[:];
σzzRes=σzzRes[isnan.(σzzRes).==0];#removing nans to get proper mean
σxzRes=σxz.-σxyAn[:];
σxzRes=σxzRes[isnan.(σxzRes).==0];#removing nans to get proper mean

using Statistics
P=[Statistics.mean(σxxRes[:]) Statistics.mean(σzzRes[:]) Statistics.mean(σxzRes[:])];
@info P
if any(P.>0.15)
 error("Error in surrounding grid is too large. This no longer matches analytical solutions, check the comparative images")
else
 println("Everything looks good, checks the mean σxx σzz σxz residual of all points is below 0.15")
end

#=
#Draw
min = minimum(σxyAn[:]);
max = maximum(σxyAn[:]);
Value=reshape(σxz,size(xobs))
x=xobs;
y=yobs;
using NaNMath
steps=10; #Steps from centre to top. 
Spacing=(max-min)./10;
levels = [min:Spacing:max;]
using PyPlot
close()
contourf(x,y,Value, levels=levels);
cbar = colorbar()
 =#