#Test case comparing to Penny shaped crack

#Start creating vars for function: 
println("creating func vars")


#Load triangles and points from file (mesh)
SurfaceDir=CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"CircleMesh_1a_500Faces.ts")
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


############Shearing############			
#Traction vector
Tn=zeros(n,1);
Tds=zeros(n,1);
Tss=ones(n,1);
FixedEls=zeros(n,1);
#Set BoundaryConditions
BoundaryConditions=Tractions(Tn,Tss,Tds)

#Calculate slip on faces
(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);

#Compute total shearing
TotalShearing = sqrt.((Dss).^2 .+(Dds).^2);


## Vol Reduction:
#Numerical
(Area,HalfPerimeter ) = CutAndDisplaceJulia.AreaOfTriangle3D( P1[:,1],P1[:,2],P1[:,3],P2[:,1],P2[:,2],P2[:,3],P3[:,1],P3[:,2],P3[:,3] );
Vol1NumResultShr=(sum(Area.*Dss)./2);
#Analytical
Radius=1;
Vol1AnResultShr=(8*(1-ν)*Tss[1]*Radius^3)/(3*(2-ν)*G);
#Equation in PPR=
Eq1=(100/Vol1AnResultShr)*(Vol1AnResultShr-Vol1NumResultShr);

println("Volume shearing - analytical then numerical then Perc difference")
println(Vol1AnResultShr)
println(Vol1NumResultShr)
println(Eq1)


############Opening############	
#Traction vector
Tn=ones(n);
Tds=zeros(n);
Tss=zeros(n);
#Set BoundaryConditions
BoundaryConditions=Tractions(Tn,Tss,Tds)

#Calculate slip on faces
(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);

## Vol Reduction (one wall):
#Numerical
Vol1NumResultOpn=(sum(Area.*Dn)./2);
#Analytical
Vol1AnResultOpn=(((pi*(1-ν)*Tn[1]*(Radius^2))/(2*G))*(2/pi))*(1/0.75)
#Vol1AnResultOpn=-(4*Tn[1]*Radius^3*(ν - 1))/(3*G);

#Equation in PPR=
Eq2=(100/Vol1AnResultOpn)*(Vol1AnResultOpn-Vol1NumResultOpn);

println("Volume Opening - analytical then numerical then Perc difference")
println(Vol1AnResultOpn)
println(Vol1NumResultOpn)
println(Eq2)

#Finding distance from 0,0 for code midpoints for plotting. 
(θ,ρ) = CutAndDisplaceJulia.cart2pol(MidPoint[:,1],MidPoint[:,2]);

#Compute analytical solution for shearing of cracks walls
a=1.0; #Radius of penny
ts=1.0;
tn=1.0
(us)=CutAndDisplaceJulia.Eshelby1957_PennyCrackSlipProfile(G,ν,0.0,ts,ρ,a)#Shearing
(un)=CutAndDisplaceJulia.Eshelby1957_PennyCrackSlipProfile(G,ν,tn,0.0,ρ,a)#Opening

#Compute the percent error between analytical and numerical
ResidualPercentDs=CutAndDisplaceJulia.BAsPercentOfA(us,TotalShearing);
ResidualPercentDn=CutAndDisplaceJulia.BAsPercentOfA(un,Dn);

#Only check values that are not within 10% of the tipline
Good=ρ.<0.9;
ResidualPercentDs=ResidualPercentDs[Good];
ResidualPercentDn=ResidualPercentDn[Good];
#@info ResidualPercentDs ResidualPercentDn

#Test this has not changed 
lim=11; #Percent error limit
if any((abs.(ResidualPercentDs.-100)).>lim) || any((abs.(ResidualPercentDn.-100)).>lim) 
 	error("Residual displacement too high, some over $lim%")
end

println("Test Passed")


#=
#To Draw
using Plots
gr()
y=[us TotalShearing];
plot1=scatter(ρ,y,title="R vs Shear Disp, An=y1 BEM=y2")
y=[un Dn];
plot2=scatter(ρ,y,title="R vs Normal Disp, An=y1 BEM=y2")
plot(plot1,plot2,layout=(2,1))
=#
using UnicodePlots
y=[us TotalShearing];
plt=scatterplot(vec(ρ),vec(us), title = "Slip of penny shaped cracks walls \n r=$a [m] ts=$ts [MPa] G=$G [MPa] ν=$ν",
 name = "analytical", xlabel = "x [m]", ylabel = "ds [m]", canvas = DotCanvas)
scatterplot!(plt, vec(ρ),vec(TotalShearing), color = :blue, name = "numerical $n tris")
println(plt) #need when running as test case

y=[us TotalShearing];
plt=scatterplot(vec(ρ),vec(un), title = "Opening of penny shaped cracks walls \n r=$a [m] tn=$tn [MPa] G=$G [MPa] ν=$ν",
 name = "analytical", xlabel = "x [m]", ylabel = "dn [m]", canvas = DotCanvas)
scatterplot!(plt, vec(ρ),vec(Dn), color = :blue, name = "numerical $n tris")
println(plt) #need when running as test case


