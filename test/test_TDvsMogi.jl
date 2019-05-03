#Test case comparing to Penny shaped crack

#Start creating vars for function: 
println("creating func vars")

SurfaceDir=CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"SphereUniformDistributionON_46Faces.ts")
(Points,Triangles)=CutAndDisplaceJulia.GoCadAsciiReader(SurfaceDir)

SurfaceDir=CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"SphereFix4.ts")
(PointsF,TrianglesF)=CutAndDisplaceJulia.GoCadAsciiReader(SurfaceDir)
PointsF[:,2:4]=PointsF[:,2:4].*3;
#Append to the two meshes
(Points,Triangles,FixedEls) = CutAndDisplaceJulia.DataAppender3D( Points,PointsF,Triangles,TrianglesF);

#Put 50m below free Surface
Depth=50;
freesurface_height=Depth;
# Getting surface ready for script. 
Points[:,4]=Points[:,4].-freesurface_height;

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
HSFlag=1; #const 

#Traction vector
Tn=ones(n);
Tds=zeros(n);
Tss=zeros(n);

#Set BoundaryConditions
BoundaryConditions=Tractions(Tn,Tss,Tds);

#Calculate slip on faces
(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);

(FixedEls,Triangles,FaceNormalVector,MidPoint,P1,P2,P3)=
CutAndDisplaceJulia.RemoveFixedElements3D(FixedEls,Triangles,FaceNormalVector,MidPoint,P1,P2,P3)

Difference=maximum(Dn)-minimum(Dn);
#@info Difference
if Difference>(0.05)/G
	error("Residual is too high")
end

P=Tn[1];
Radius=1;
X =[0.0:1:100.0;]
Y=zeros(size(X));
Z=zeros(size(X));
X=reshape(X,length(X),1);
Y=reshape(Y,length(Y),1);
Z=reshape(Z,length(Z),1);



#Compute analytical result
(UxAn,UyAn,UzAn)=CutAndDisplaceJulia.Mogi1958_SphericalCavity(Depth,X,Y,Radius,P,ν,G)

DispFlag=1;
StrainFlag=0;

StrainInfVector=Strains([],[],[],[],[],[]);
DispInfVector=Disps([],[],[]);
(StrainAtInfPoints,DispAtInfPoints)= 
CutAndDisplaceJulia.TD(X,Y,Z,P1,P2,P3,Dss,Dds,Dn,ν,G,DispFlag,StrainFlag,HSFlag,StrainInfVector,DispInfVector)


#Compute the percent error between analytical and numerical
ResidualPercentUx=CutAndDisplaceJulia.BAsPercentOfA(UxAn,DispAtInfPoints.Ux);
ResidualPercentUy=CutAndDisplaceJulia.BAsPercentOfA(UyAn,DispAtInfPoints.Uy);
ResidualPercentUz=CutAndDisplaceJulia.BAsPercentOfA(UzAn,DispAtInfPoints.Uz);

#@info ResidualPercentUx
#@info ResidualPercentUy
#@info ResidualPercentUz

lim=20; #Percent error limit
if any((abs.(ResidualPercentUx.-100)).>lim) | any((abs.(ResidualPercentUy.-100)).>lim)  | any((abs.(ResidualPercentUz.-100)).>lim) 
 	error("Residual displacement too high, some over $lim%")
end

println("Test Passed")


#To Draw
#using Plots
#gr()
#y=[UxAn.+UyAn.+UzAn Ux.+Uy.+Uz];
#scatter(X,y,title="X vs TotalDisp, An=y1 BEM=y2")
using UnicodePlots
y=[UxAn.+UyAn.+UzAn DispAtInfPoints.Ux.+DispAtInfPoints.Uy.+DispAtInfPoints.Uz];
plt=lineplot(vec(X),vec(y[:,1]), title = "Deformation above mogi source \n r=$Radius [m] P=$P [MPa] G=$G [MPa] ν=$ν", name = "analytical", xlabel = "x [m]", ylabel = "ux+uy+uz [m]", canvas = DotCanvas)
lineplot!(plt, vec(X),vec(y[:,2]), color = :blue, name = "numerical $n tris")
println(plt) #need when running as test case
