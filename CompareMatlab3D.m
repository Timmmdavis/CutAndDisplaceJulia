%Comment - start some vectors (spaced points)

[X,Y,Z] = meshgrid(-3:.02:3,-3:.002:3,2);
X=1;
Y=1;
Z=1;

P1=[-1 0 0];
P2=[1 -1 -1];
P3=[0 1.5 .5];

Ss=-1;
Ds=2;
Ts=3;
nu=0.25;

tic
[U]=TDdispFS(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,nu);
toc
%disp(U)