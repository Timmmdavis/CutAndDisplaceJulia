clear

%Fault geom
Strike=60;
Dip=45;
Length=5;
Width=5;
TipDepth=1;
Rake=90;
%Fault slip
Dds=2;
Dn=0.6;
Dss=0.5;
nu=0.25;
mu=1;
[ ~,E,lambda,nu,mu ] = ElasticConstantsCheck( mu,nu );

%%%Section - arranging the fault surface for the TDE solution
%Arranged as so
%P1=(-1  1 0 ; 1 -1 0); %top left and bottom right
%P2=(-1 -1 0 ;-1 -1 0); %bottom left
%P3=( 1  1 0 ; 1  1 0); %top right
X=[-1  1 -1 -1 1 1];
Y=[ 1 -1 -1 -1 1 1];
Z=[ 0  0  0  0 0 0];
X=(X./2).*Width;
Y=(Y./2).*Length;

%Getting the direction cosines for the actual plane
[ StrikeSlipCosine,DipSlipCosine ] = CalculateDSandSSDirs( [cosd(Dip),0,cosd(Dip)] ) ;%Plane dipping at dip (facing east!)
RAngle=deg2rad(-Strike);
StrikeSlipCosine=RotateCosine3d(StrikeSlipCosine,RAngle,"z");
DipSlipCosine=RotateCosine3d(DipSlipCosine,RAngle,"z");
FaceNormalVector=cross((StrikeSlipCosine),(DipSlipCosine));
FaceNormalVector=FaceNormalVector';


%Now rotate the flat plane to the correct location. 
[X,Y,Z]=RotateObject3dNewCoords(DipSlipCosine,StrikeSlipCosine,FaceNormalVector,X(:),Y(:),Z(:));
Drop=max(Z);
Z=Z-Drop-TipDepth; %Move fault down
% % P1=[X(1) Y(1) Z(1);X(2) Y(2) Z(2)];
% % P2=[X(5) Y(5) Z(5);X(4) Y(4) Z(4)];
% % P3=[X(3) Y(3) Z(3);X(6) Y(6) Z(6)];

P1=[X(1) Y(1) Z(1);X(2) Y(2) Z(2)];
P2=[X(5) Y(5) Z(5);X(4) Y(4) Z(4)];
P3=[X(3) Y(3) Z(3);X(6) Y(6) Z(6)];

%disp('increasing no of tris')
%P1=repmat(P1,50,1);
%P2=repmat(P2,50,1);
%P3=repmat(P3,50,1);

% Start some vectors (spaced points)
x = linspace(-10,10,50);
[x,y]=meshgrid(x,x);
z=ones(size(x))*-0; %Ground surface

%Get lengths (for reshapes later)
[dimx,dimy] = size(x);
%Turn to col vectors
x=reshape(x,numel(x),1);
y=reshape(y,numel(y),1);
z=reshape(z,numel(z),1);

halfspace=1;
FD=0;

X=x(:);
Y=y(:);
Z=z(:);

tic
NUM=size(X,1);
Stressinfmatrix = zeros(NUM,6); 
Dispinfmatrix=[];
Dss=1;  
Dds=0;  
Dn=0; 
[Dssinfmatrix,DssDisplacementXYZ]=CreateCoeffsLoop3d(Stressinfmatrix,Dispinfmatrix,...
NUM,X,Y,Z,P1,P2,P3,Dss,Dds,Dn,mu,lambda,nu,halfspace,FD);
toc



% x=reshape(x,dimx,dimy);
% y=reshape(y,dimx,dimy);
% z=reshape(z,dimx,dimy);
% ux=reshape(DisplacementXYZ(:,1),dimx,dimy);
% contourf(x,y,ux);