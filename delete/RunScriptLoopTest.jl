#Start creating vars for function: 
println("creating func vars")

# Comment - start some vectors (spaced points)
x = [-2:0.05:2;];
y = [-4:0.05:0;];
x,y=CutAndDisplaceJulia.meshgrid(x,y);
#Get lengths (for reshapes later)
dimx,dimy = size(x);
#Turn to col vectors
x=reshape(x,length(x),1);
y=reshape(y,length(y),1);

#Other vars 
xe=0;
ye=-2;
a=1a;
Beta=deg2rad(0);

R1= Array{Float64}(undef, 1,length(x));
R2= Array{Float64}(undef, 1,length(x));

#Just calling by self
@time (R1,R2)=CutAndDisplaceJulia.DotLoopTestVect(x,y,xe,ye,a,Beta);
@time (R1,R2)=CutAndDisplaceJulia.DotLoopTest(x,y,xe,ye,a,Beta);