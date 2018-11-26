#Test case comparing Okada dislocations with TDs (displacement)
#Would be better to manually rotate points etc (not just predefine)

#Start creating vars for function: 
println("creating func vars")

#Inputs

#Fault geom
Strike=60;
Dip=70;
Length=5;
Width=5;
TipDepth=1;
Rake=90;
#Fault slip
Dds=1;
Dn=0;
Dss=0;
#Elastic cons
nu=0.25;

#Arranged as so
#P1=[-1  1 0 ; 1 -1 0]; #top left and bottom right
#P2=[-1 -1 0 ;-1 -1 0]; #bottom left
#P3=[ 1  1 0 ; 1  1 0]; #top right

X=[-1  1 -1 -1 1 1];
Y=[ 1 -1 -1 -1 1 1];
Z=[ 0  0  0  0 0 0];

X=X.*Width;
Y=Y.*Length;

( StrikeSlipCosine,DipSlipCosine ) = MyModule.CalculateDSandSSDirs( cosd(45),0,cosd(45) ) #Plane dipping at dip (facing east!)
#THEN ROTATE IT BY STRIKE (NEED COSINE ROTATE METHO)
(X,Y,Z)=MyModule.RotateObject3DNewCoords(X,Y,Z,0,0,0,Ax1,Ax2,Ax3)

#Define dipping fault plane (Would be better to add rectangle and transform....
#P1=[2.59258869000000	0.509504670000000	-5.69846310000000; -2.59258869000000	-0.509504670000000	-1]
#P2=[1.73753833000000	1.99049533000000	-1; -1.73753833000000	-1.99049533000000	-5.69846310000000]
#P3=[-1.73753833000000	-1.99049533000000	-5.69846310000000; 1.73753833000000	1.99049533000000	-1]


# Comment - start some vectors (spaced points)
x = range(-10,stop=10,length=50); #linspace deprecated
(x,y)=MyModule.meshgrid(x,x);
z=ones(size(x))*-0; #Ground surface

#Get lengths (for reshapes later)
dimx,dimy = size(x);
#Turn to col vectors
x=reshape(x,length(x),1);
y=reshape(y,length(y),1);
z=reshape(z,length(z),1);




println("Vars created -> to TD func1")
(Ux1,Uy1,Uz1)=MyModule.TDdispHS(x,y,z,P1[1,:],P2[1,:],P3[1,:],Dss,Dds,Dn,nu);
(Ux2,Uy2,Uz2)=MyModule.TDdispHS(x,y,z,P1[2,:],P2[2,:],P3[2,:],Dss,Dds,Dn,nu);
Ux=Ux1.+Ux2;
Uy=Uy1.+Uy2;
Uz=Uz1.+Uz2;
println("Out of func, too Okada")

MidDeptha=sind(Dip)*Width;
MidDepth=MidDeptha/2+TipDepth;
println("Into Okada func")
(uX,uY,uZ)=MyModule.Okada1985RectangularDislocation(x,y,MidDepth,Strike,Dip,Length,Width,Rake,Dds,Dn,nu);

println("Out of func, drawing time, start by reshape")
x=reshape(x,dimx,dimy);
y=reshape(y,dimx,dimy);
uX=reshape(uX,dimx,dimy);
Ux=reshape(Ux,dimx,dimy);

UxRes=maximum(Ux[:].-uX[:]);
UyRes=maximum(Uy[:].-uY[:]);
UzRes=maximum(Uz[:].-uZ[:]);


if UxRes>1E-9
	error("UxRes too high, Okada and TD not matching for displacement")
end
if UyRes>1E-9
	error("UyRes too high, Okada and TD not matching for displacement")
end
if UzRes>1E-9
	error("UzRes too high, Okada and TD not matching for displacement")
end

println("Test Passed")

#= if you want to draw remove this line
println("Starting Makie")
using Makie

println("Calling Func 2 draw")
sc1=contour(x, y, uX, levels = 20, linewidth = 0, fillrange = true)
#sc2=contour(x, y, Ux, levels = 20, linewidth = 0, fillrange = true)
#hbox(sc2, sc1)
=#