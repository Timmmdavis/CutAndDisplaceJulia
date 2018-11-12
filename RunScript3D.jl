#Start creating vars for function: 
println("creating func vars")

# Comment - start some vectors (spaced points)
X = [-3:0.02:3;];
Y = [-3:0.002:3;];
X,Y=MyModule.meshgrid(X,Y);
Z=ones(size(X))*2;

#Get lengths (for reshapes later)
dimx,dimy = size(X);
#Turn to col vectors
X=reshape(X,length(X),1);
Y=reshape(Y,length(Y),1);

P1=[-1,0,0];
P2=[1,-1,-1];
P3=[0,1.5,.5];

Ss=-1;
Ds=2;
Ts=3;
nu=0.25;

X=1;
Y=1;
Z=1;

println("Vars created -> to func")

tic=time()
(ue,un,uv)=MyModule.TDdispFS(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,nu);
toc=time()
println("Elapsed time")
println(toc-tic)

#using Profile
#@profile  (ue,un,uv)=MyModule.TDdispFS(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,nu);
##Profile.print()
#Profile.print(format=:flat)