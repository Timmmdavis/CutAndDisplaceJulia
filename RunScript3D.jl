#Start creating vars for function: 
println("creating func vars")

# Comment - start some vectors (spaced points)
x = [-3:0.02:3;];
y = [-3:0.02:3;];
x,y=MyModule.meshgrid(x,y);
z=ones(size(x))*2;
#Get lengths (for reshapes later)
dimx,dimy = size(x);
#Turn to col vectors
x=reshape(x,1,length(x));
y=reshape(y,1,length(y));
z=reshape(z,1,length(z));
#Points
P1=[-1.,0.,0.];
P2=[1.,-1.,-1.];
P3=[0.,1.5,.5];
#Angles
Ax=pi;
Ay=0;
Az=0;

Ss=-1;
Ds=2;
Ts=3;
nu=0.25;
mu=6;
lambda=(2*mu*nu)/(1-(2*nu));

#X=1;
#Y=1;
#Z=1;

println("Vars created -> to func")

ue=zeros(size(x));
un=zeros(size(x));
uv=zeros(size(x));


#Timing
#@time (ue,un,uv)=MyModule.TDdispFS(x,y,z,P1,P2,P3,Ss,Ds,Ts,nu);
@time (Exx,Eyy,Ezz,Exy,Exz,Eyz)=MyModule.TDstrainFS(x,y,z,P1,P2,P3,Ss,Ds,Ts,mu,lambda);
un=Exx;

#using Profile
#@profile  (ue,un,uv)=MyModule.TDdispFS(x,y,z,P1,P2,P3,Ss,Ds,Ts,nu);
#Profile.print(format=:filefuncline )
#Profile.print(format=:flat)
#Profile.print(sortedby=:count)

#using MyModule
#using Traceur
#@trace MyModule.TDdispFS(x,y,z,P1,P2,P3,Ss,Ds,Ts,nu);

#=
tic=time()
Threads.@threads for i=1:100 #
	(ue[:,1],un[:,1],uv[:,1])=MyModule.TDdispFS(x,y,z,P1,P2,P3,Ss,Ds,Ts,nu);
	#println(i)
end
toc=time()
println("Elapsed time")
println(toc-tic)
=#

####Some reshaping for drawing:
#Was one func that reshaped stuff
x=reshape(x,dimx,dimy);
y=reshape(y,dimx,dimy);
un=reshape(un,dimx,dimy);

#Draw
using NaNMath
Top=maximum([NaNMath.maximum(un),abs(NaNMath.minimum(un))])
steps=10; #Steps from centre to top. 
levels = [-Top:Top/steps:Top;]
using PyPlot
close()
contourf(x,y,un, levels=levels);
cbar = colorbar()
