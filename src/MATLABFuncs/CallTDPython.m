%Create Arrays and convert to NumPy
clear

%% Define MATLAB matricies
x = linspace(-10,10,4);
[x,y]=meshgrid(x,x);
z=ones(size(x))*-0; %Ground surface
x=reshape(x,numel(x),1);
y=reshape(y,numel(y),1);
z=reshape(z,numel(z),1);

X=[-1  1 -1 -1 1 1];
Y=[ 1 -1 -1 -1 1 1];
Z=[ 0  -3  0  -3 0 0];
P1=[X(1) Y(1) Z(1);X(2) Y(2) Z(2)];
P2=[X(5) Y(5) Z(5);X(4) Y(4) Z(4)];
P3=[X(3) Y(3) Z(3);X(6) Y(6) Z(6)];

DdsVec=[2.0,2.0];
DnVec=[0.6,0.6];
DssVec=[0.5,0.5];

mu=1.; %const 
nu=0.25; %const 

DispFlag=1;
StressFlag=1;
HSFlag=0;

%% Convert to Python 
[x,y,z] = ConvertMatriciesToNumPy(x,y,z);

[P1,P2,P3] = ConvertMatriciesToNumPy(P1,P2,P3);

[DdsVec,DnVec,DssVec] = ConvertMatriciesToNumPy(DdsVec,DnVec,DssVec);

mu=py.float(mu);
nu=py.float(nu);

DispFlag=py.int(DispFlag);
StressFlag=py.int(StressFlag);
HSFlag=py.int(HSFlag); 

%% Call script
py.RunJuliaTD.TDInterface(x,y,z,P1,P2,P3,DssVec,DdsVec,DnVec,nu,mu,DispFlag,StressFlag,HSFlag);

%% Clear and load results
clear
load ('PythonOutput.mat')