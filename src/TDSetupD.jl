function TDSetupD(x,y,z,alpha,bx,by,bz,nu,TriVertex,SideVec)
# TDSetupD transforms coordinates of the calculation points as well as 
# slip vector components from ADCS into TDCS. It then calculates the 
# displacements in ADCS and transforms them into TDCS.

(Ct,St,y1,z1,by1,bz1)=TransformToADCS(y,z,by,bz,SideVec,TriVertex)

#Init arrays
Ang=-pi+alpha;
cosA = cos(Ang);
sinA = sin(Ang);
uV  = Array{Float64}(undef, length(x),1);
v0V = Array{Float64}(undef, length(x),1);
w0V = Array{Float64}(undef, length(x),1);
#Extra defs out of loop to speed it up
E1=(1-nu); #Elastic cons
E2=(1-2*nu);
cosA2=cosA^2;
sinADE1=sinA/8/pi/(1-nu);
# Calculate displacements associated with an angular dislocation in ADCS
for i=1:length(x)

	(u,v0,w0) = AngDisDisp(x[i],y1[i],z1[i],cosA,sinA,bx,by1[1], bz1[1], E1,E2,cosA2,sinADE1);

	uV[i]=u;
	v0V[i]=v0;
	w0V[i]=w0;
end

# Transform displacements from ADCS into TDCS
(v,w)  =RotateObject2D(v0V,w0V,0,0,Ct,-St) #Rotate back

return(uV,v,w)

end