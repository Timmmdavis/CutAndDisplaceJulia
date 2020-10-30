#Test case comparing to Penny shaped crack

#Elastic constants
G=ShearModulus(2.0e9); 
ν=PoissonsRatio(0.25);#println("Pr is close to 0.5")
(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(G,ν);

#Density
ρrock=2900;
ρfluid=2600;
#Gravitational acceleration
g=9.81;

Δρ=((ρrock-ρfluid)*g)

KCrit=5e7; #[5e7 = 50 MPa √m]

# Top tip critically stressed
CritRadius=abs(real((complex(9*Δρ)^(1/3)*(pi*(KCrit^2))^(1/3)))/(4*Δρ))

P=(2*Δρ*CritRadius)/3


#Load triangles and points from file (mesh)
SurfaceDir=CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"CircleMesh_1a_500Faces.off")
(Points,Triangles)=CutAndDisplaceJulia.OFFReader(SurfaceDir)

#=
#Flatten such that normals point up
X=Points[:,2];
Points[:,2]=Points[:,4];
Points[:,4]=X;
Radius=maximum(Points[:,2])
=#

#Semiaxes= (a must be bigger than b (tada solution))
a=CritRadius; 
b=CritRadius;
Points[:,2]=Points[:,2]*a;
Points[:,3]=Points[:,3]*b;
#Put vertical
#Pennys angle away from Z. 
Beta=0; #45 
BetaFromVert=90-Beta;
#Rotate this (YZ)
(Points[:,3],Points[:,4])=CutAndDisplaceJulia.RotateObject2D!(Points[:,3],Points[:,4],0.0,0.0,cosd(BetaFromVert),sind(BetaFromVert))



#Define a number of tris you want themesh to have
(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )
(target_edge_length,max_target_edge_length)=
CutAndDisplaceJulia.GetDesiredEdgeLength(P1,P2,P3,1500) #650

#Remesh using Polygon method in CGAL:
OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"BeforePolygonRemshing")
(OutputDirectory)=BuildCGAL.PolygonRemeshingCGAL(OutputDirectory,target_edge_length)
(Points,Triangles)=CutAndDisplaceJulia.OFFReader(OutputDirectory)
(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.CleanEdgeTris(Points,Triangles)
(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.IsosceliseEdgeTrisNew(MidPoint,P1,P2,P3,Triangles,Points,FaceNormalVector)

P1,P2,P3,MidPoint,FaceNormalVector,Points,Triangles=
CutAndDisplaceJulia.PutOnEdgeOfCirc(MidPoint,P1,P2,P3,FaceNormalVector,CritRadius,1,3)

#=
using Plots
pyplot()
plot4 = plot()
P1nrm=P1./a;
P2nrm=P2./a;
P3nrm=P3./a;
for i=1:length(P1nrm[:,1])
    Plots.plot!([P1nrm[i,1],P2nrm[i,1]],[P1nrm[i,3],P2nrm[i,3]], aspect_ratio=:equal,c=(:black), lab="")
	Plots.plot!([P2nrm[i,1],P3nrm[i,1]],[P2nrm[i,3],P3nrm[i,3]], aspect_ratio=:equal,c=(:black), lab="")
	Plots.plot!([P1nrm[i,1],P3nrm[i,1]],[P1nrm[i,3],P3nrm[i,3]], aspect_ratio=:equal,c=(:black), lab="")
end 
display(plot4)
=#

#Do Twice:
println("Doing twice")
OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"BeforePolygonRemshing")
(OutputDirectory)=BuildCGAL.PolygonRemeshingCGAL(OutputDirectory,target_edge_length)
(Points,Triangles)=CutAndDisplaceJulia.OFFReader(OutputDirectory)
(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.CleanEdgeTris(Points,Triangles)
(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.IsosceliseEdgeTrisNew(MidPoint,P1,P2,P3,Triangles,Points,FaceNormalVector)

P1,P2,P3,MidPoint,FaceNormalVector,Points,Triangles=
CutAndDisplaceJulia.PutOnEdgeOfCirc(MidPoint,P1,P2,P3,FaceNormalVector,CritRadius,1,3)

println("Doing 3x")
OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"BeforePolygonRemshing")
(OutputDirectory)=BuildCGAL.PolygonRemeshingCGAL(OutputDirectory,target_edge_length)
(Points,Triangles)=CutAndDisplaceJulia.OFFReader(OutputDirectory)
P1,P2,P3,MidPoint,FaceNormalVector,Points,Triangles=
CutAndDisplaceJulia.PutOnEdgeOfCirc(MidPoint,P1,P2,P3,FaceNormalVector,CritRadius,1,3)
(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.CleanEdgeTris(Points,Triangles)
(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.IsosceliseEdgeTrisNew(MidPoint,P1,P2,P3,Triangles,Points,FaceNormalVector)


OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"OutMesh")

n=length(Triangles[:,1]);
n2=length(Points[:,1]);

# Which bits we want to compute
HSFlag=0; #const 

Z=MidPoint[:,3];

#Set BoundaryConditions
σxx = zeros(n);  
σyy = zeros(n);  
σzz = zeros(n);  
σxy = zeros(n);    
σxz = zeros(n);
σyz = zeros(n);
#Bouyancy
Tn=zeros(n)
Tss=zeros(n)
Tds=zeros(n)

for j=1:length(Z)
	Tn[j]=Δρ.*Z[j]; 
	Tn[j]=Tn[j]+P;
end


Stress=Stresses(σxx,σyy,σzz,σxy,σxz,σyz);
Traction=Tractions(Tn,Tss,Tds);
BoundaryConditions=MixedBoundaryConditions(Stress,Traction)

FixedEls=zeros(n,1);

#Calculate slip on faces
(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);

#Get tip elements
(FeP1P2S,FeP1P3S,FeP2P3S)=CutAndDisplaceJulia.GetCrackTipElements3D(MidPoint,P1,P2,P3,FaceNormalVector)
(FeP1P2S,FeP1P3S,FeP2P3S)=CutAndDisplaceJulia.StressIntensity3D(Dn,Dss,Dds,G,ν,FaceNormalVector,FeP1P2S,FeP1P3S,FeP2P3S);

#Calculate theta [location around the crack
FreeEdMdX=[FeP1P2S.FeMd[FeP1P2S.FreeFlg,1];FeP1P3S.FeMd[FeP1P3S.FreeFlg,1];FeP2P3S.FeMd[FeP2P3S.FreeFlg,1]];
FreeEdMdY=[FeP1P2S.FeMd[FeP1P2S.FreeFlg,2];FeP1P3S.FeMd[FeP1P3S.FreeFlg,2];FeP2P3S.FeMd[FeP2P3S.FreeFlg,2]];
FreeEdMdZ=[FeP1P2S.FeMd[FeP1P2S.FreeFlg,3];FeP1P3S.FeMd[FeP1P3S.FreeFlg,3];FeP2P3S.FeMd[FeP2P3S.FreeFlg,3]];

θnum=atan.(FreeEdMdX,FreeEdMdZ);

#Compute across each theta
K1an,θ=CutAndDisplaceJulia.Tada_StrIntEllipseCrackTension(P,a,b,FreeEdMdX,FreeEdMdY)
K1anx=CutAndDisplaceJulia.Tada_StrIntPennyGradient(Δρ,a,θnum)
K1an=K1an.+K1anx; #sum these


avgeverynth=3
FeP1P2S,FeP1P3S,FeP2P3S=CutAndDisplaceJulia.MovingAverageOfStressIntensity(avgeverynth,P1,P2,P3,FaceNormalVector,MidPoint,FeP1P2S,FeP1P3S,FeP2P3S)

#%Accumulate from structure into big vectors:
K1=[FeP1P2S.K1[FeP1P2S.FreeFlg];FeP1P3S.K1[FeP1P3S.FreeFlg];FeP2P3S.K1[FeP2P3S.FreeFlg]];
K2=[FeP1P2S.K2[FeP1P2S.FreeFlg];FeP1P3S.K2[FeP1P3S.FreeFlg];FeP2P3S.K2[FeP2P3S.FreeFlg]];
K3=[FeP1P2S.K3[FeP1P2S.FreeFlg];FeP1P3S.K3[FeP1P3S.FreeFlg];FeP2P3S.K3[FeP2P3S.FreeFlg]];


#Compute the percent error between analytical and numerical
ResidualPercentK1=CutAndDisplaceJulia.BAsPercentOfA(K1an,K1);

X=K1an-K1
println("vertical sep")
aaa=maximum(abs.(K1an.-K1))./maximum(K1an)
println(aaa)
L2Norm=sqrt(sum(abs.(X).^2))
println("L2norm")
@info L2Norm

#To Draw
IntAng=[FeP1P2S.IntAng[FeP1P2S.FreeFlg];FeP1P3S.IntAng[FeP1P3S.FreeFlg];FeP2P3S.IntAng[FeP2P3S.FreeFlg]];
#mutable struct TriangleEdges;FeLe;FeMd;FeEv;FeM2Ev;FreeFlg;FeM2ELe;IntAng;K1;K2;K3;StrainEnergy; 	end
Area=[FeP1P2S.Area[FeP1P2S.FreeFlg];FeP1P3S.Area[FeP1P3S.FreeFlg];FeP2P3S.Area[FeP2P3S.FreeFlg]];


Area=Area./IntAng

using Plots
pyplot()


#plot1=scatter(θ,K1an,c=(:black))
#using LaTeXStrings
t=-180:1:180;
x = a.*cosd.(t) 
y = a.*sind.(t) 
K1any,zz=CutAndDisplaceJulia.Tada_StrIntEllipseCrackTension(P,a,b,x,y)
K1anx=CutAndDisplaceJulia.Tada_StrIntPennyGradient(Δρ,a,deg2rad.(t))
K1an=K1any.+K1anx; #sum these
plot3=Plots.plot(t,K1an./maximum(K1an),c=(:black), lab=latexstring(raw"$K_{I}$"))
scatter!(rad2deg.(θnum),K1./maximum(K1an),c=(:black),ms=1.5, lab=latexstring("BEM"))
yaxis!("Stress intensity", (0,1.5))
xaxis!(latexstring("\$\\theta^{\\circ}\$ angle around crack"), (-180,180))
Plots.title!(latexstring("Penny crack, \$K_{Ic}=1\$"))
#plot1=scatter(θ,K1an,'k')
#scatter!(θ,K1,zcolor=ResidualPercentK1, m=(:reds),markersize=(Area./maximum(Area)).*10)
#plot2=scatter(FreeEdMdX,FreeEdMdY,markersize=(Area./maximum(Area)).*10, aspect_ratio=:equal,zcolor=ResidualPercentK1, m=(:reds), lab="")

plot4 = Plots.plot()

P1nrm=P1./a;
P2nrm=P2./a;
P3nrm=P3./a;

for i=1:length(P1[:,1])

    Plots.plot!([P1nrm[i,1],P2nrm[i,1]],[P1nrm[i,3],P2nrm[i,3]], aspect_ratio=:equal,c=(:black), lab="")
	Plots.plot!([P2nrm[i,1],P3nrm[i,1]],[P2nrm[i,3],P3nrm[i,3]], aspect_ratio=:equal,c=(:black), lab="")
	Plots.plot!([P1nrm[i,1],P3nrm[i,1]],[P1nrm[i,3],P3nrm[i,3]], aspect_ratio=:equal,c=(:black), lab="")
end 

yaxis!(latexstring("\$z\$"), (-1.5,1.5))
xaxis!(latexstring("\$x\$"), (-1.5,1.5))
annotate!(0, 1.1, text(latexstring("\$\\theta\$=0\$^{\\circ}\$"), 10, :center))
annotate!(-1.4, 0, text(latexstring("\$\\theta\$=-90\$^{\\circ}\$"), 10, :left))
annotate!(1.4,  0, text(latexstring("\$\\theta\$=90\$^{\\circ}\$"), 10, :right))


PLT=plot(plot3,plot4,layout=(2,1))
display(PLT)

#If you want all plots
#PLT=plot(plot2,plot1,plot4,plot3,layout=(2,2))

#=
#@info ResidualPercentK1 ResidualPercentK2 ResidualPercentK3
MaxErrorK1=maximum(filter(!isnan,abs.(ResidualPercentK1))); 
@info MaxErrorK1 
#Test this has not changed 
lim=11; #Percent error limit
if MaxErrorK1>lim
 	error("Residual displacement too high, some over $lim%")
end
println("Test Passed")
=#



#=
σzz=σzz[1]
using UnicodePlots
plt=scatterplot(vec(θ),vec(K1an), title = "StressIntensity around inclined penny subject to tension \n r=$Radius [m] σzz=$σzz [MPa] βFromZ=$BetaFromVert [°] G=$G [MPa] ν=$ν", name = "analytical KI", xlabel = "θ [°]", ylabel = "K", canvas = DotCanvas)
scatterplot!(plt, vec(θ),vec(K1), color = :blue, name = "numerical KI $n tris")
println(plt) #need when running as test case
plt=scatterplot(vec(θ),vec(K2an), title = "StressIntensity around inclined penny subject to tension \n r=$Radius [m] σzz=$σzz [MPa] βFromZ=$BetaFromVert [°]  G=$G [MPa] ν=$ν", name = "analytical KII", xlabel = "θ [°]", ylabel = "K", canvas = DotCanvas)
scatterplot!(plt, vec(θ),K2, color = :yellow, name = "numerical KII $n tris")
println(plt) #need when running as test case
plt=scatterplot(vec(θ),vec(K3an), title = "StressIntensity around inclined penny subject to tension \n r=$Radius [m] σzz=$σzz [MPa] βFromZ=$BetaFromVert [°]  G=$G [MPa] ν=$ν", name = "analytical KIII", xlabel = "θ [°]", ylabel = "K", canvas = DotCanvas)
scatterplot!(plt, vec(θ),vec(K3), color = :magenta, name = "numerical KIII $n tris")
println(plt) #need when running as test case
=#

