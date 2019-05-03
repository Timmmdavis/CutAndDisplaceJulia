
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

############################################################
println("Surface normal along y")
#Now make surface vertical (Switch Y and Z) strikes at 90
PointsTmp=Points[:,4];
Points[:,4]=Points[:,3];
Points[:,3]=PointsTmp;
#Compute triangle properties
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
n=length(Triangles[:,1]);
n2=length(Points[:,1]);
(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )

println("Surface normal along x")
#Now make surface vertical (Switch X and Y) strikes at 0
PointsTmp=Points[:,3];
Points[:,3]=Points[:,2];
Points[:,2]=PointsTmp;

#Roatating so no dip.
#=
Ang=1; #degree
(Points[:,2],Points[:,4])=CutAndDisplaceJulia.RotateObject2D!( Points[:,2],Points[:,4],0.,0.,cosd(Ang),sind(Ang) )
=#

#Compute triangle properties
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
n=length(Triangles[:,1]);
n2=length(Points[:,1]);

(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )

CutAndDisplaceJulia.CleanAndIsosceliseEdgeTris(MidPoint,P1,P2,P3,Triangles,FaceNormalVector)

############################################################

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

BoundaryConditions=Stresses(σxx,σyy,σzz,ones(n),σxz,σyz);
#Calculate slip on faces

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
Ang=[Ang1; Ang2; Ang3]

a=2; #plot for x
b=3; #plot for z
XMid=[FeP1P2S.FeMd[FeP1P2S.FreeFlg,a]
	  FeP1P3S.FeMd[FeP1P3S.FreeFlg,a]
	  FeP2P3S.FeMd[FeP2P3S.FreeFlg,a]]

YMid=[FeP1P2S.FeMd[FeP1P2S.FreeFlg,b]
	  FeP1P3S.FeMd[FeP1P3S.FreeFlg,b]
	  FeP2P3S.FeMd[FeP2P3S.FreeFlg,b]]

XDir=[FeP1P2S.FeEv[FeP1P2S.FreeFlg,a]
	  FeP1P3S.FeEv[FeP1P3S.FreeFlg,a]
	  FeP2P3S.FeEv[FeP2P3S.FreeFlg,a]]

YDir=[FeP1P2S.FeEv[FeP1P2S.FreeFlg,b]
	  FeP1P3S.FeEv[FeP1P3S.FreeFlg,b]
	  FeP2P3S.FeEv[FeP2P3S.FreeFlg,b]]

XDir2=[FeP1P2S.FeM2Ev[FeP1P2S.FreeFlg,a]
	   FeP1P3S.FeM2Ev[FeP1P3S.FreeFlg,a]
	   FeP2P3S.FeM2Ev[FeP2P3S.FreeFlg,a]]

YDir2=[FeP1P2S.FeM2Ev[FeP1P2S.FreeFlg,b]
	   FeP1P3S.FeM2Ev[FeP1P3S.FreeFlg,b]
	   FeP2P3S.FeM2Ev[FeP2P3S.FreeFlg,b]]	

K2=[FeP1P2S.K2[FeP1P2S.FreeFlg,1]
    FeP1P3S.K2[FeP1P3S.FreeFlg,1]
    FeP2P3S.K2[FeP2P3S.FreeFlg,1]]	

K1=[FeP1P2S.K1[FeP1P2S.FreeFlg,1]
    FeP1P3S.K1[FeP1P3S.FreeFlg,1]
    FeP2P3S.K1[FeP2P3S.FreeFlg,1]]	

K3=[FeP1P2S.K3[FeP1P2S.FreeFlg,1]
    FeP1P3S.K3[FeP1P3S.FreeFlg,1]
    FeP2P3S.K3[FeP2P3S.FreeFlg,1]]  

strainenergy=[FeP1P2S.StrainEnergy[FeP1P2S.FreeFlg,1]
              FeP1P3S.StrainEnergy[FeP1P3S.FreeFlg,1]
              FeP2P3S.StrainEnergy[FeP2P3S.FreeFlg,1]]          	     


#Functionize and return final fig (all 3) - or just group before?
using Plots
gr()
scl=0.1;#length of vectors
fig = plot(reuse=false,legend=false)
for i=1:length(XMid)

	x=XMid[i]; 
	y=YMid[i]; 
	scatter!([x],[y],ms=strainenergy[i]*1e5) #ms=abs(K3[i])*4
	xnew=x+XDir[i]*scl; 
	ynew=y+YDir[i]*scl; 
	xnew2=x+XDir2[i]*scl; 
	ynew2=y+YDir2[i]*scl; 
	xd1=[x xnew]
	yd1=[y ynew]
	xd2=[x xnew2]
	yd2=[y ynew2]
	plot!(vec(xd1),vec(yd1))
	plot!(vec(xd2),vec(yd2))

end 	



println("σxy∞ twist test")
XDir=[FeP1P2S.FeM2Ev[FeP1P2S.FreeFlg,1]
	  FeP1P3S.FeM2Ev[FeP1P3S.FreeFlg,1]
	  FeP2P3S.FeM2Ev[FeP2P3S.FreeFlg,1]]
YDir=[FeP1P2S.FeM2Ev[FeP1P2S.FreeFlg,2]
	  FeP1P3S.FeM2Ev[FeP1P3S.FreeFlg,2]
	  FeP2P3S.FeM2Ev[FeP2P3S.FreeFlg,2]]
ZDir=[FeP1P2S.FeM2Ev[FeP1P2S.FreeFlg,3]
	  FeP1P3S.FeM2Ev[FeP1P3S.FreeFlg,3]
	  FeP2P3S.FeM2Ev[FeP2P3S.FreeFlg,3]]
#2nd we check edges pointing along neg y direction
flag=YDir.==-1
#because the top tips of this edge dont progate in the expected direction of a 2D fracture:
#We compute a ratio between those propagating properly and those not. 
ratio=sum(Px[flag].>0)/sum(Px[flag].<0)
println("Good and bad point twist:")
@info sum(Px[flag].>0) sum(Px[flag].<0)
using Plots
gr()
scl=0.1;#length of vectors
fig2 = plot(reuse=false,legend=false)
scatter!(Px,Pz,ms=-Py*4) #Looking towards pos y
@info Px

GoodEdges=findall(Px[flag].>0);
BadEdges=findall(Px[flag].<0);
@info GoodEdges BadEdges K2

#Get the index's of the points:
Freelocs=findall([FeP1P2S.FreeFlg; FeP1P3S.FreeFlg; FeP2P3S.FreeFlg])
@info size(Freelocs)
Indxs=findall(flag)

#A list of index's to the edges that are pointing towards the south
#if over 'n' then the next set. i.e. i:n = P1P2 n+1:2n=P1P3 2n+1:3n=P2P3
Indxes=Freelocs[Indxs]
@info Indxes Indxes[BadEdges] n

GoodIndx=Indxes[GoodEdges];
BadIndx=Indxes[BadEdges];

#Find index locations:
IndxP1P2=findall(FeP1P2S.FreeFlg);
IndxP1P3=findall(FeP1P3S.FreeFlg);
IndxP2P3=findall(FeP2P3S.FreeFlg);
@info IndxP1P2

Fe=FeP1P2S;
for i=1:length(IndxP1P2)
    #Extract some values
    I=IndxP1P2[i]; 
    if any(I.==BadIndx)
    	println("Bad == P1P2 at $i")
    	println(FaceNormalVector[I,:])
    	println(Fe.FeMd[I,:])
    	println(Fe.FeM2Ev[I,:])
        println(Fe.FeLe[I])
        println(Fe.FeEv[I,:])
        println(Fe.K2[I])
        println(Fe.K1[I])
    end
end

Fe=FeP1P3S;
for i=1:length(IndxP1P3)
    #Extract some values
    I=IndxP1P3[i]; 
    if any(I.==(BadIndx.-n))
    	println("Bad == P1P3 at $i")
    	println(FaceNormalVector[I,:])
    	println(Fe.FeMd[I,:])
    	println(Fe.FeM2Ev[I,:])
        println(Fe.FeLe[I])
        println(Fe.FeEv[I,:])
        println(Fe.K2[I])
        println(Fe.K1[I])
    end
end

Fe=FeP2P3S;
for i=1:length(IndxP2P3)
    #Extract some values
    I=IndxP2P3[i]; 
    if any(I.==(BadIndx.-2*n))
    	println("Bad == P2P3 at $i")
    	println(FaceNormalVector[I,:])
    	println(Fe.FeMd[I,:])
    	println(Fe.FeM2Ev[I,:])
        println(Fe.FeLe[I])
        println(Fe.FeEv[I,:])
        println(Fe.K2[I])
        println(Fe.K1[I])
    end
end


K2_InterestedEdge=K2[flag];
K1_InterestedEdge=K1[flag];

#Part of vectors
XDir=[FeP1P2S.FeEv[FeP1P2S.FreeFlg,1]
	  FeP1P3S.FeEv[FeP1P3S.FreeFlg,1]
	  FeP2P3S.FeEv[FeP2P3S.FreeFlg,1]]
YDir=[FeP1P2S.FeEv[FeP1P2S.FreeFlg,2]
	  FeP1P3S.FeEv[FeP1P3S.FreeFlg,2]
	  FeP2P3S.FeEv[FeP2P3S.FreeFlg,2]]
ZDir=[FeP1P2S.FeEv[FeP1P2S.FreeFlg,3]
	  FeP1P3S.FeEv[FeP1P3S.FreeFlg,3]
	  FeP2P3S.FeEv[FeP2P3S.FreeFlg,3]]
XDir_InterestedEdge=XDir[flag];
YDir_InterestedEdge=YDir[flag];
ZDir_InterestedEdge=ZDir[flag];

Ang_InterestedEdge=Ang[flag];

@info K2_InterestedEdge[GoodEdges] K1_InterestedEdge[GoodEdges]
@info K2_InterestedEdge[BadEdges]  K1_InterestedEdge[BadEdges]

@info XDir_InterestedEdge[GoodEdges] YDir_InterestedEdge[GoodEdges] ZDir_InterestedEdge[GoodEdges]
@info XDir_InterestedEdge[BadEdges]  YDir_InterestedEdge[BadEdges] ZDir_InterestedEdge[BadEdges]

@info Ang[GoodEdges]  
@info Ang[BadEdges]  

return fig



