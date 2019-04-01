#Test case comparing to 2D crack with linearly increasing cohesive profile.

#Start creating vars for function: 
println("creating func vars")

#Load triangles and points from file (mesh)
SurfaceDir=CutAndDisplaceJulia.LoadData(CutAndDisplaceJulia,"BladeFaultShort_500Tris_Vertical_NorthDip_0.5.ts")
(Points,Triangles)=CutAndDisplaceJulia.GoCadAsciiReader(SurfaceDir)

println("Flattening surface in Y")
Points[:,3]=Points[:,3]*0; 
println("Scaling on")
SclCrck=0.0001;
Points[:,2:4]=Points[:,2:4].*SclCrck;

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
σxx = 0.0;	
σyy = -0.05*SclCrck;#force faces closed	
σzz = 0.0;	
σxy = 1*SclCrck;	
σxz = 0.0;
σyz = 0.0;

###########WithFriction###################
Stress=Stresses(σxx,σyy,σzz,σxy,σxz,σyz);
FixedEls=zeros(n,1);
#Set BoundaryConditions
Traction=Tractions(0.0,0.0,0.0)
#Friction parameters
µ=0.0;
Sf  = abs.(MidPoint[:,1]); 
Friction=FrictionParameters(µ,Sf);
BoundaryConditions=MixedBoundaryConditionsFriction(Stress,Traction,Friction);
#Calculate slip on faces
(DnF, DssF, DdsF)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);

###########WithoutFriction###################
BoundaryConditions=Stresses(σxx,σyy,σzz,σxy,σxz,σyz);
#Calculate slip on faces
(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);


#Driving Traction
Ts=abs.(σxy);
a = SclCrck;  #Unit half length 
x=MidPoint[:,1]; #Along x axis
Good=abs.(MidPoint[:,3]).<(SclCrck); #Remove points that are not close to crack centre
x=x[Good]
Dss=abs.(Dss[Good])
DssF=abs.(DssF[Good])
#Analytical solution (no Friction)
(Slip_UniformRemote,~)=CutAndDisplaceJulia.PollardSegall1987_FractureSlipProfile(G,ν,0.0,Ts,x,a);
#Analytical solution (Friction)
Slip_IncreasingFriction=CutAndDisplaceJulia.Burgmann1994_FractureLinearFrictionSlipProfile(G,ν,Ts,x,a);

#Compute the percent error between analytical and numerical
ResidualPercentDs=CutAndDisplaceJulia.BAsPercentOfA(Slip_UniformRemote,Dss);
ResidualPercentDsF=CutAndDisplaceJulia.BAsPercentOfA(Slip_IncreasingFriction,DssF);

ResidualSumDs=sum(ResidualPercentDs.-100)
ResidualSumDsF=sum(ResidualPercentDsF.-100)
@info ResidualSumDs ResidualSumDsF

#Only check values that are not within 10% of the tipline
Good=(abs.(x)).<0.9*SclCrck;
ResidualPercentDs=ResidualPercentDs[Good];
ResidualPercentDsF=ResidualPercentDsF[Good];
@info ResidualPercentDs ResidualPercentDsF

#Test this has not changed 
lim=36; #Percent error limit
if any((abs.(ResidualPercentDs.-100)).>lim) || any((abs.(ResidualPercentDsF.-100)).>lim) 
 	error("Residual displacement too high, some over $lim%")
end
lim=1400; #Sum of percent error
if ResidualSumDs>lim || ResidualSumDsF>lim
 	error("Sum of Residual displacement percentages too high, over $lim%")
end


#=
#To Draw
using Plots
gr()
y=[Slip_UniformRemote Dss];
plot1=scatter(x,y,title="R vs Shear Disp, An=y1 BEM=y2")
y=[Slip_IncreasingFriction DssF];
plot2=scatter(x,y,title="R vs Normal Disp, An=y1 BEM=y2")
plot(plot1,plot2,layout=(2,1))
=#


###################################################################################################
#Run flat surface

PointsTmp=Points[:,4];
Points[:,4]=Points[:,3];
Points[:,3]=PointsTmp;
println("Scaling on")
SclCrck2=1e4;
Points[:,2:4]=Points[:,2:4]./SclCrck; #reset to 1 as half length
SclCrck=SclCrck2;
Points[:,2:4]=Points[:,2:4].*SclCrck; #use a difference scaling

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
σxx = 0.0;	
σyy = 0.0;#force faces closed	
σzz = -0.05*SclCrck;	
σxy = 0.0;	
σxz = 1*SclCrck;
σyz = 0.0;

###########WithFriction###################
Stress=Stresses(σxx,σyy,σzz,σxy,σxz,σyz);
FixedEls=zeros(n,1);
#Set BoundaryConditions
Traction=Tractions(0.0,0.0,0.0)
#Friction parameters
µ=0.0;
Sf  = abs.(MidPoint[:,1]); 
Friction=FrictionParameters(µ,Sf);
BoundaryConditions=MixedBoundaryConditionsFriction(Stress,Traction,Friction);
#Calculate slip on faces
(DnF, DssF, DdsF)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);

###########WithoutFriction###################
BoundaryConditions=Stresses(σxx,σyy,σzz,σxy,σxz,σyz);
#Calculate slip on faces
(Dn, Dss, Dds)=CutAndDisplaceJulia.SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls);


#Driving Traction
Ts=abs.(σxz);
a = SclCrck;  #Unit half length 
x=MidPoint[:,1]; #Along x axis
Good=abs.(MidPoint[:,2]).<(SclCrck); #Remove points that are not close to crack centre
x=x[Good]
Dds=abs.(Dds[Good])
DdsF=abs.(DdsF[Good])
#Analytical solution (no Friction)
(Slip_UniformRemote,~)=CutAndDisplaceJulia.PollardSegall1987_FractureSlipProfile(G,ν,0.0,Ts,x,a);
#Analytical solution (Friction)
Slip_IncreasingFriction=CutAndDisplaceJulia.Burgmann1994_FractureLinearFrictionSlipProfile(G,ν,Ts,x,a);

#Compute the percent error between analytical and numerical
ResidualPercentDs=CutAndDisplaceJulia.BAsPercentOfA(Slip_UniformRemote,Dds);
ResidualPercentDsF=CutAndDisplaceJulia.BAsPercentOfA(Slip_IncreasingFriction,DdsF);

ResidualSumDs=sum(ResidualPercentDs.-100)
ResidualSumDsF=sum(ResidualPercentDsF.-100)
@info ResidualSumDs ResidualSumDsF

#Only check values that are not within 10% of the tipline
Good=(abs.(x)).<0.9*SclCrck;
ResidualPercentDs=ResidualPercentDs[Good];
ResidualPercentDsF=ResidualPercentDsF[Good];
@info ResidualPercentDs ResidualPercentDsF

#Test this has not changed 
lim=36; #Percent error limit
if any((abs.(ResidualPercentDs.-100)).>lim) || any((abs.(ResidualPercentDsF.-100)).>lim) 
 	error("Residual displacement too high, some over $lim%")
end
lim=1400; #Sum of percent error
if ResidualSumDs>lim || ResidualSumDsF>lim
 	error("Sum of Residual displacement percentages too high, over $lim%")
end


#=
#To Draw
using Plots
gr()
y=[Slip_UniformRemote Dds];
plot1=scatter(x,y,title="R vs Shear Disp, An=y1 BEM=y2")
y=[Slip_IncreasingFriction DdsF];
plot2=scatter(x,y,title="R vs Normal Disp, An=y1 BEM=y2")
plot(plot1,plot2,layout=(2,1))
=#


println("Test Passed")