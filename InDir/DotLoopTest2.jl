function DotLoopTest2(X,Y,Z,P1,P2,P3,Ss,Ds,Ts,nu)
# TDdispFS 

#=
bx = Ts; # Tensile-slip
by = Ss; # Strike-slip
bz = Ds; # Dip-slip

#Init some vars
p1 = zeros(3,1);
p2 = zeros(3,1);
p3 = zeros(3,1);

eY = [0;1;0];
eZ = [0;0;1];
P1_1=P1[1];P1_2=P1[2];P1_3=P1[3];
P2_1=P2[1];P2_2=P2[2];P2_3=P2[3];
P3_1=P3[1];P3_2=P3[2];P3_3=P3[3];
=#

# Calculate unit strike, dip and normal to TD vectors: For a horizontal TD 
# as an exception, if the normal vector points upward, the strike and dip 
# vectors point Northward and Westward, whereas if the normal vector points
# downward, the strike and dip vectors point Southward and Westward, 
# respectively.
P2P1= Array{Float64}(undef, size(P1));
P3P1= Array{Float64}(undef, size(P1));
for i=1:length(P1)
	P2P1[i] = P2[i]-P1[i];
	P3P1[i] = P2[i]-P1[i];
end
Vnorm = cross(P2P1,P3P1);
NVnorm=norm(Vnorm);
Vstrike = cross(eZ,Vnorm);
NVstrike=norm(Vstrike)
if NVstrike==0
	for i=1:length(Vstrike)
		Vstrike = eY[i]*Vnorm[3];
	end
end
NVstrike=norm(Vstrike)
for i=1:length(Vnorm)
	Vnorm[i] = Vnorm[i]/NVnorm[i];
	Vstrike[i] = Vstrike[i]/NVstrike[i];
end
Vdip = cross(Vnorm,Vstrike);

#uelp = Array{Float64}(undef, 1,length(X)); 
#unlp = Array{Float64}(undef, 1,length(X)); 
#uvlp = Array{Float64}(undef, 1,length(X)); 

end #remove

