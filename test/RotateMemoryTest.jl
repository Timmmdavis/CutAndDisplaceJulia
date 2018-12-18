########################################################## Rotate

# # Start some vectors (spaced points)
# x = range(-10,stop=10,length=5000); #linspace deprecated
# (x,y)=MyModule.meshgrid(x,x);
# z=ones(size(x))*-0; #Ground surface

# #Get lengths (for reshapes later)
# dimx,dimy = size(x);
# #Turn to col vectors
# x=reshape(x,length(x),1);
# y=reshape(y,length(y),1);
# z=reshape(z,length(z),1);

# th=deg2rad(45)
# Ct=cos(th);
# St=sin(th);

# Ax1=[1 0 0]
# Ax2=[0 -1 0]
# Ax3=[0 0 1]

# MyModule.TestFooAllocs(x,y,z,0,0,0,Ct,St,Ax1,Ax2,Ax3)


########################################################## Tens trans

theta=45.;
Ct=cosd(theta);
St=sind(theta);

# #println("Remove Allocation here")
# A = [[1 0  0  ];
     # [0 Ct St ];
	 # [0 -St Ct]]; # 3x3 Transformation matrix
# B = [[1 0  0  ];
     # [0 Ct -St ];
	 # [0 St Ct]]; # 3x3 Transformation matrix


Ax1=[1   0   0 ]
Ax2=[0  Ct  St ]
Ax3=[0  -St  Ct]	 
	 
A=[Ax1; Ax2; Ax3]	 
B=[[Ax1[1] Ax2[1] Ax3[1]]; [Ax1[2] Ax2[2] Ax3[2]]; [Ax1[3] Ax2[3] Ax3[3]]]
	 
Exx=rand(100000,1);
Eyy=rand(100000,1);
Ezz=rand(100000,1);
Exy=rand(100000,1);
Exz=rand(100000,1);
Eyz=rand(100000,1);

MyModule.TestFooAllocs2(Exx,Eyy,Ezz,Exy,Exz,Eyz,A,B)
