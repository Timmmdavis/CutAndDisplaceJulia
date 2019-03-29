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
Tn=1.0;
Tds=0.0;
Tss=0.0;

#Set BoundaryConditions
BoundaryConditions=Tractions(Tn,Tss,Tds);

#Calculate slip on faces
(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);

println("PutInFunc")
#Get logical flag of bits to keep
FD=FixedEls.==0;
#Needs to be a single vector
FD=reshape(FD,length(FD))
#Remove bits. 
Triangles=Triangles[FD,:];
FaceNormalVector=FaceNormalVector[FD,:];
MidPoint=MidPoint[FD,:];
P1=P1[FD,:];
P2=P2[FD,:];
P3=P3[FD,:];



Difference=maximum(Dn)-minimum(Dn);
@info Difference
if Difference>(0.05)/G
	error("Residual is too high")
end

P=Tn;
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
StressFlag=0;
#Compute displacements
(ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn,
 ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
 ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,
 UxDn,UyDn,UzDn,
 UxDss,UyDss,UzDss,
 UxDds,UyDds,UzDds)=
 CutAndDisplaceJulia.TD(X,Y,Z,P1,P2,P3,Dss,Dds,Dn,ν,G,DispFlag,StressFlag,HSFlag)
 

(Exx,Eyy,Ezz,Exy,Exz,Eyz,Ux,Uy,Uz)=
CutAndDisplaceJulia.TD_sum(ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
	   ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,
	   UxDn,UyDn,UzDn,
	   UxDss,UyDss,UzDss,
       UxDds,UyDds,UzDds)


#Compute the percent error between analytical and numerical
ResidualPercentUx=CutAndDisplaceJulia.BAsPercentOfA(UxAn,Ux);
ResidualPercentUy=CutAndDisplaceJulia.BAsPercentOfA(UyAn,Uy);
ResidualPercentUz=CutAndDisplaceJulia.BAsPercentOfA(UzAn,Uz);

@info ResidualPercentUx
@info ResidualPercentUy
@info ResidualPercentUz

lim=20; #Percent error limit
if any((abs.(ResidualPercentUx.-100)).>lim) | any((abs.(ResidualPercentUy.-100)).>lim)  | any((abs.(ResidualPercentUz.-100)).>lim) 
 	error("Residual displacement too high, some over $lim%")
end

println("Test Passed")


#To Draw
# using Plots
# gr()
# y=[UxAn.+UyAn.+UzAn Ux.+Uy.+Uz];
# scatter(X,y,title="X vs TotalDisp, An=y1 BEM=y2")