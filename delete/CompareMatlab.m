%CompareMatlab

[x,y]=meshgrid(-2:0.05:2,-4:0.05:0);
dimx = length(x(:,1));
dimy = length(x(1,:));
Ds=0; Dn=1;

x=x(:);
y=y(:);
xe=0;
ye=-2;
a=1;
Beta=0;
nu=0.25;
E=1;


Stress=zeros(size(x));
Stress=[Stress,Stress,Stress];

% tic
% for i=1
%     for j=1:numel(x)
%         [Stress] = LDstressHS(x(j),y(j),xe,ye,a,Beta,Ds,Dn,nu,E);
%     end
%     disp(i)
% end
% toc


tic
for i=1:500
    [Stress] = LDstressHS(x,y,xe,ye,a,Beta,Ds,Dn,nu,E);
    disp(i)
end
toc