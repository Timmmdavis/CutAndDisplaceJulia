#Start creating vars for function: 
println("creating func vars")

# Comment - start some vectors (spaced points)
x = [-2:0.005:2;];
y = [-4:0.05:0;];
x,y=MyModule.meshgrid(x,y);
#Get lengths (for reshapes later)
dimx,dimy = size(x);
#Turn to col vectors
x=reshape(x,length(x),1);
y=reshape(y,length(y),1);

#Other vars 
Ds=0;
Dn=1;
xe=0;
ye=-2;
a=1;
Beta=deg2rad(0);
nu=0.25;
E=1;

#Just calling by self
#(Sxx,Syy,Sxy)=MyModule.LDstressHS(x,y,xe,ye,a,Beta,Ds,Dn,nu,E);
#println(Sxx)

println("Vars created -> to func")
#Init array (stress at each point) (Only needed for loops). 
Sxx= Array{Float64}(undef, 1,length(x));
Syy= Array{Float64}(undef, 1,length(x));
Sxy= Array{Float64}(undef, 1,length(x));

####Profiling
#=
###Run Func Pre profile
x2=[2;2;2;2]; y2=[2;2;2;2]; 
Sxx2= Array{Float64}(undef, 1,length(x));Syy2=Sxx2;Szz2=Sxx2;
(Sxx2,Syy2,Sxy2)=MyModule.LDstressHS(x2,y2,xe,ye,a,Beta,Ds,Dn,nu,E);

using Profile
 #Threads.@threads
 @profile Sxx[1,:],Syy[1,:],Sxy[1,:]=MyModule.LDstressHS(x,y,xe,ye,a,Beta,Ds,Dn,nu,E)
Profile.print()
Profile.print(format=:flat)


####Start loop
tic=time()
Threads.@threads for i=1:500 #

	(Sxx[1,:],Syy[1,:],Sxy[1,:])=MyModule.LDstressHS(x,y,xe,ye,a,Beta,Ds,Dn,nu,E);

	#@show i
end
toc=time()
println("Elapsed time")
println(toc-tic)


println("Out of function")

####Some reshaping for drawing:
#Was one func that reshaped stuff
x=reshape(x,dimx,dimy);
y=reshape(y,dimx,dimy);
Sxx=reshape(Sxx,dimx,dimy);
Syy=reshape(Syy,dimx,dimy);
Sxy=reshape(Sxy,dimx,dimy);

print("Off to PyPlot")

#Draw
using NaNMath
Top=maximum([NaNMath.maximum(Sxx),abs(NaNMath.minimum(Sxx))])
steps=10; #Steps from centre to top. 
levels = [-Top:Top/steps:Top;]
using PyPlot
contourf(x, y,Sxx, levels=levels);
cbar = colorbar()
=#