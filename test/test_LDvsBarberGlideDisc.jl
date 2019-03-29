#Start creating vars for function: 
println("creating func vars")

#Inputs
G=ShearModulus(500.); 
ν=PoissonsRatio(0.25);
(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(G,ν);


#create x and y
spacing=0.1;
minx=-4; maxx=4;
x = [minx:spacing:maxx;];
(X,Y) = CutAndDisplaceJulia.meshgrid(x,x);
#Get lengths (for reshapes later)
dimx,dimy = size(X);

#Doing with two Els
a=1;
hlflen = [0.5 0.5];  
xe=[-0.5 0.5];
ye=[0. 0.];
Beta=[0. 0.];
Ds=[1. 1.];
Dn=[0. 0.];


DispFlag=1;
StressFlag=1;
HSFlag=0;

println("Vars created -> to LD func")



(σxxDs,σyyDs,σxyDs,σxxDn,σyyDn,σxyDn,UxDs,UxDn,UyDs,UyDn)=CutAndDisplaceJulia.LD(X[:],Y[:],xe,ye,hlflen,Beta,Ds,Dn,ν,G,DispFlag,StressFlag,HSFlag)
#Accumulating arrays
(σxx,σyy,σxy,Ux,Uy)=CutAndDisplaceJulia.LD_sum(σxxDs,σxxDn,σyyDs,σyyDn,σxyDs,σxyDn,UxDs,UxDn,UyDs,UyDn)

println("Out of func, too Barber1992_GlideDislocation")
(X,Y,σxxAn,σyyAn,σxyAn,UxAn,UyAn)=CutAndDisplaceJulia.Barber1992_GlideDislocation(G,X[:],Y[:],a,Ds[1],ν)


println("Out of func, drawing time, start by reshape")

UxRes=maximum(filter(!isnan, UxAn[:].-Ux[:]));
UyRes=maximum(filter(!isnan, UyAn[:].-Uy[:]));
σxxRes=maximum(filter(!isnan, σxxAn[:].-σxx[:]));
σyyRes=maximum(filter(!isnan, σyyAn[:].-σyy[:]));
σxyRes=maximum(filter(!isnan, σxyAn[:].-σxy[:]));

println("Values of residuals: LD vs Barber")
@info UxRes UyRes σxxRes σyyRes σxyRes
if UxRes>1E-13
	error("UxRes too high, LD and glide dislocation not matching for displacement")
end
if UyRes>1E-13
	error("UyRes too high, LD and glide dislocation not matching for displacement")
end
if σxxRes>1E-12
	error("SxxRes too high, LD and glide dislocation not matching for stress")
end
if σyyRes>1E-12
	error("SyyRes too high, LD and glide dislocation not matching for stress")
end
if σxyRes>1E-12
	error("SxyRes too high, LD and glide dislocation not matching for stress")
end

println("Test Passed")

#=
#Draw
struct PlotAn;	σxx;σyy;σxy; end
struct PlotNum;	σxx;σyy;σxy; end
AnalyticalResults=PlotAn(σxxAn,σyyAn,σxyAn)
NumericalResults=PlotNum(σxx,σzz,σxz)
#Get fields in structure
FieldsInStructAn=fieldnames(typeof(AnalyticalResults));
FieldsInStructNum=fieldnames(typeof(NumericalResults));
for i=1:length(FieldsInStructAn)
	AnStr=FieldsInStructAn[i];
	NumStr=FieldsInStructNum[i];
	println(AnStr)
	AnVal=getfield(AnalyticalResults, FieldsInStructAn[i])
	NumVal=getfield(NumericalResults, FieldsInStructNum[i])
	min = minimum(AnVal[:]);
	max = maximum(AnVal[:]);
	NumValue=reshape(NumVal,size(xobs))
	AnValue=reshape(AnVal,size(xobs))
	x=xobs;
	y=yobs;
	using NaNMath
	steps=10; #Steps from centre to top. 
	Spacing=(max-min)./10;
	levels = [min:Spacing:max;]
	using PyPlot
	figure()
	plt.subplot(211)
	contourf(x,y,NumValue, levels=levels);
	cbar = colorbar();
	plt.title("Numerical Result $AnStr");
	plt.subplot(212)
	contourf(x,y,AnValue, levels=levels);
	cbar = colorbar();
	plt.title("Analytical Result $NumStr")
end
 =#

