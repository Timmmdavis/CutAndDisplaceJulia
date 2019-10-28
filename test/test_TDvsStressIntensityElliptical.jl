#Test case comparing to Penny shaped crack

#Start creating vars for function: 
println("creating func vars")

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
a=2; 
b=1;
Points[:,2]=Points[:,2]*a;
Points[:,3]=Points[:,3]*b;
CritRadius=1

#Define a number of tris you want the mesh to have
(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )#650
(target_edge_length,max_target_edge_length)=
CutAndDisplaceJulia.GetDesiredEdgeLength(P1,P2,P3,1500)#650

#Remesh using Polygon method in CGAL:
OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"BeforePolygonRemshing")
(OutputDirectory)=BuildCGAL.PolygonRemeshingCGAL(OutputDirectory,target_edge_length)
(Points,Triangles)=CutAndDisplaceJulia.OFFReader(OutputDirectory)
(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.CleanEdgeTris(Points,Triangles)
(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.IsosceliseEdgeTrisNew(MidPoint,P1,P2,P3,Triangles,Points,FaceNormalVector)


#Put back into circle
Points[:,2]=Points[:,2]./a;
Points[:,3]=Points[:,3]./b;
(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
P1,P2,P3,MidPoint,FaceNormalVector,Points,Triangles=
CutAndDisplaceJulia.PutOnEdgeOfCirc(MidPoint,P1,P2,P3,FaceNormalVector,CritRadius,1,2)
#Back 2 ellipse
Points[:,2]=Points[:,2].*a;
Points[:,3]=Points[:,3].*b;
(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)


#Do Twice:
println("Doing twice")
OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"BeforePolygonRemshing")
(OutputDirectory)=BuildCGAL.PolygonRemeshingCGAL(OutputDirectory,target_edge_length)
(Points,Triangles)=CutAndDisplaceJulia.OFFReader(OutputDirectory)
(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.CleanEdgeTris(Points,Triangles)
(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.IsosceliseEdgeTrisNew(MidPoint,P1,P2,P3,Triangles,Points,FaceNormalVector)


#Put back into circle
Points[:,2]=Points[:,2]./a;
Points[:,3]=Points[:,3]./b;
(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)
P1,P2,P3,MidPoint,FaceNormalVector,Points,Triangles=
CutAndDisplaceJulia.PutOnEdgeOfCirc(MidPoint,P1,P2,P3,FaceNormalVector,CritRadius,1,2)
#Back 2 ellipse
Points[:,2]=Points[:,2].*a;
Points[:,3]=Points[:,3].*b;
(P1,P2,P3)=CutAndDisplaceJulia.CreateP1P2P3( Triangles,Points )
(FaceNormalVector,MidPoint)=CutAndDisplaceJulia.CreateFaceNormalAndMidPoint(Points,Triangles)

println("Doing 3x")
OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"BeforePolygonRemshing")
(OutputDirectory)=BuildCGAL.PolygonRemeshingCGAL(OutputDirectory,target_edge_length)
(Points,Triangles)=CutAndDisplaceJulia.OFFReader(OutputDirectory)
(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.CleanEdgeTris(Points,Triangles)
(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.IsosceliseEdgeTrisNew(MidPoint,P1,P2,P3,Triangles,Points,FaceNormalVector)
#=
println("Doing 4x")
OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"BeforePolygonRemshing")
(OutputDirectory)=BuildCGAL.PolygonRemeshingCGAL(OutputDirectory,target_edge_length)
(Points,Triangles)=CutAndDisplaceJulia.OFFReader(OutputDirectory)
(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.CleanEdgeTris(Points,Triangles)
(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.IsosceliseEdgeTrisNew(MidPoint,P1,P2,P3,Triangles,Points,FaceNormalVector)
println("Doing 5x")
OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"BeforePolygonRemshing")
(OutputDirectory)=BuildCGAL.PolygonRemeshingCGAL(OutputDirectory,target_edge_length)
(Points,Triangles)=CutAndDisplaceJulia.OFFReader(OutputDirectory)
(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.CleanEdgeTris(Points,Triangles)
(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.IsosceliseEdgeTrisNew(MidPoint,P1,P2,P3,Triangles,Points,FaceNormalVector)
println("Doing 6x")
OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"BeforePolygonRemshing")
(OutputDirectory)=BuildCGAL.PolygonRemeshingCGAL(OutputDirectory,target_edge_length)
(Points,Triangles)=CutAndDisplaceJulia.OFFReader(OutputDirectory)
(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.CleanEdgeTris(Points,Triangles)
(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.IsosceliseEdgeTrisNew(MidPoint,P1,P2,P3,Triangles,Points,FaceNormalVector)
println("Doing 7x")
OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"BeforePolygonRemshing")
(OutputDirectory)=BuildCGAL.PolygonRemeshingCGAL(OutputDirectory,target_edge_length)
(Points,Triangles)=CutAndDisplaceJulia.OFFReader(OutputDirectory)
(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.CleanEdgeTris(Points,Triangles)
(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.IsosceliseEdgeTrisNew(MidPoint,P1,P2,P3,Triangles,Points,FaceNormalVector)
println("Doing 8x")
OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"BeforePolygonRemshing")
(OutputDirectory)=BuildCGAL.PolygonRemeshingCGAL(OutputDirectory,target_edge_length)
(Points,Triangles)=CutAndDisplaceJulia.OFFReader(OutputDirectory)
(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.CleanEdgeTris(Points,Triangles)
(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.IsosceliseEdgeTrisNew(MidPoint,P1,P2,P3,Triangles,Points,FaceNormalVector)
println("Doing 9x")
OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"BeforePolygonRemshing")
(OutputDirectory)=BuildCGAL.PolygonRemeshingCGAL(OutputDirectory,target_edge_length)
(Points,Triangles)=CutAndDisplaceJulia.OFFReader(OutputDirectory)
(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.CleanEdgeTris(Points,Triangles)
(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.IsosceliseEdgeTrisNew(MidPoint,P1,P2,P3,Triangles,Points,FaceNormalVector)
println("Doing 10x")
OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"BeforePolygonRemshing")
(OutputDirectory)=BuildCGAL.PolygonRemeshingCGAL(OutputDirectory,target_edge_length)
(Points,Triangles)=CutAndDisplaceJulia.OFFReader(OutputDirectory)
(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.CleanEdgeTris(Points,Triangles)
(P1,P2,P3,Triangles,Points,MidPoint,FaceNormalVector)=CutAndDisplaceJulia.IsosceliseEdgeTrisNew(MidPoint,P1,P2,P3,Triangles,Points,FaceNormalVector)
=#


OutputDirectory=CutAndDisplaceJulia.OFFExport(Points,Triangles,length(Triangles[:,1]),length(Points[:,1]),"OutMesh")

n=length(Triangles[:,1]);
n2=length(Points[:,1]);

# Which bits we want to compute
HSFlag=0; #const 

#Elastic constants
G=ShearModulus(5.); 
ν=PoissonsRatio(0.25);
(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(G,ν);

#Set BoundaryConditions
σxx = zeros(n);	
σzz = ones(n);	
σyy = zeros(n);	
σxy = zeros(n);    
σxz = zeros(n);
σyz = zeros(n);
BoundaryConditions=Stresses(σxx,σyy,σzz,σxy,σxz,σyz);

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

#Compute across each theta
K1an,θ=CutAndDisplaceJulia.Tada_StrIntEllipseCrackTension(σzz[1],a,b,FreeEdMdX,FreeEdMdY)

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
t=[181:0.01:360; 1:0.01:180];
x = a.*cosd.(t) 
y = a.*sind.(t) 
K1an,t=CutAndDisplaceJulia.Tada_StrIntEllipseCrackTension(σzz[1],a,b,x,y)
plot1=plot(rad2deg.(t),K1an./maximum(K1an),c=(:black), lab=latexstring(raw"$K_{I}$"))
scatter!(rad2deg.(θ),K1./maximum(K1an),c=(:black),ms=1.5, lab=latexstring("BEM"))
yaxis!("Stress intensity", (0,1.5))
xaxis!(latexstring("\$\\theta^{\\circ}\$ angle around crack"), (-180,180))
title!(latexstring("Elliptical crack in tension, \$a_x/a_y\$ = $a/$b"))
#plot1=scatter(θ,K1an,'k')
#scatter!(θ,K1,zcolor=ResidualPercentK1, m=(:reds),markersize=(Area./maximum(Area)).*10)
#plot2=scatter(FreeEdMdX,FreeEdMdY,markersize=(Area./maximum(Area)).*10, aspect_ratio=:equal,zcolor=ResidualPercentK1, m=(:reds), lab="")


plot2 = plot()
for i=1:length(P1[:,1])

    Plots.plot!([P1[i,1],P2[i,1]],[P1[i,2],P2[i,2]], aspect_ratio=:equal,c=(:black), lab="")
	Plots.plot!([P2[i,1],P3[i,1]],[P2[i,2],P3[i,2]], aspect_ratio=:equal,c=(:black), lab="")
	Plots.plot!([P1[i,1],P3[i,1]],[P1[i,2],P3[i,2]], aspect_ratio=:equal,c=(:black), lab="")
end 

yaxis!(latexstring("\$y\$"), (-1.2,1.2))
xaxis!(latexstring("\$x\$"), (-2.5,2.5))
annotate!(0, 1.1, text(latexstring("\$\\theta\$=90\$^{\\circ}\$"), 10, :center))
annotate!(0, -1.1, text(latexstring("\$\\theta\$=-90\$^{\\circ}\$"), 10, :center))
annotate!(2.4, 0, text(latexstring("\$\\theta\$=0\$^{\\circ}\$"), 10, :right))


PLT=plot(plot1,plot2,layout=(2,1))
display(PLT)

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