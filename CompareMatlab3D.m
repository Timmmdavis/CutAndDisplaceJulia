%Comment - start some vectors (spaced points)

[X,Y,Z] = meshgrid(-3:.02:3,-3:.02:3,2);
[dimx,dimy] = size(X);
X=X(:);Y=Y(:);Z=Z(:);

%X=1;
%Y=1;
%Z=1;

P1=[-1 0 0];
P2=[1 -1 -1];
P3=[0 1.5 .5];

Ss=-1;
Ds=2;
Ts=3;
nu=0.25;
mu=6;
lambda=(2*mu*nu)/(1-(2*nu));

U=zeros(size(X,1),3);
D=zeros(size(X,1),6);

tic
for i=1:2
    i
    %U(:,:)=TDdispFS(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,nu);
    D(:,:)=TDstrainFS(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,mu,lambda);
end
toc
%disp(U)

[X,Y,Un ]=ReshapeData2d( dimx,dimy,X,Y,D(:,1) );
contourf(X,Y,Un)