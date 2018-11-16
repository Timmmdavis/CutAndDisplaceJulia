#Start creating vars for function: 
println("creating func vars")

# Comment - start some vectors (spaced points)
X = [-3:0.02:3;];
Y = [-3:0.02:3;];
X,Y=MyModule.meshgrid(X,Y);
Z=ones(size(X))*2;
#Get lengths (for reshapes later)
dimx,dimy = size(X);
#Turn to col vectors
X=reshape(X,length(X),1);
Y=reshape(Y,length(Y),1);
Z=reshape(Z,length(Z),1);
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

#X=1;
#Y=1;
#Z=1;

println("Vars created -> to func")

ue=zeros(size(X));
un=zeros(size(X));
uv=zeros(size(X));


#Timing
# @time (ue,un,uv)=MyModule.TDdispFSLooped(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,nu);

#using Profile
#@profile  (ue,un,uv)=MyModule.TDdispFSLooped(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,nu);
#Profile.print(format=:filefuncline )
#Profile.print(format=:flat)
#Profile.print(sortedby=:count)

#using MyModule
#using Traceur
#@trace MyModule.TDdispFSLooped(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,nu);

tic=time()
Threads.@threads for i=1:100 #
	(ue,un,uv)=MyModule.TDdispFSLooped(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,nu);
	#println(i)
end
toc=time()
println("Elapsed time")
println(toc-tic)