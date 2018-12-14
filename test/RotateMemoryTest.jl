
# Start some vectors (spaced points)
x = range(-10,stop=10,length=5000); #linspace deprecated
(x,y)=MyModule.meshgrid(x,x);
z=ones(size(x))*-0; #Ground surface

#Get lengths (for reshapes later)
dimx,dimy = size(x);
#Turn to col vectors
x=reshape(x,length(x),1);
y=reshape(y,length(y),1);
z=reshape(z,length(z),1);

th=deg2rad(45)
Ct=cos(th);
St=sin(th);

Ax1=[1 0 0]
Ax2=[0 -1 0]
Ax3=[0 0 1]

MyModule.TestFooAllocs(x,y,z,0,0,0,Ct,St,Ax1,Ax2,Ax3)