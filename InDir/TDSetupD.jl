function TDSetupD(x,y,z,alpha,bx,by,bz,nu,TriVertex,SideVec,casenLog)
# TDSetupD transforms coordinates of the calculation points as well as 
# slip vector components from ADCS into TDCS. It then calculates the 
# displacements in ADCS and transforms them into TDCS.

Ct=SideVec[3];
St=SideVec[2];
P1=TriVertex[2];
P2=TriVertex[3];
# Transform coordinates of the calculation points from TDCS into ADCS
(y1,z1)  =RotateObject2D(y,z,P1,P2,Ct,St)
# Transform the in-plane slip vector components from TDCS into ADCS
(by1,bz1)=RotateObject2D(by,bz,0,0,Ct,St)


# Configuration II (rotate other way)
if any(casenLog)
	(by1n,bz1n)=RotateObject2D(by,bz,0,0,Ct,St)
end


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
	if casenLog[i]==1
		(u,v0,w0) = AngDisDisp(x[i],y1[i],z1[i],cosA,sinA,bx,by1n[1],bz1n[1],E1,E2,cosA2,sinADE1);
	else
		(u,v0,w0) = AngDisDisp(x[i],y1[i],z1[i],cosA,sinA,bx,by1[1],bz1[1],E1,E2,cosA2,sinADE1);
	end
	uV[i]=u;
	v0V[i]=v0;
	w0V[i]=w0;
end

# Transform displacements from ADCS into TDCS
(v,w)  =RotateObject2D(v0V,w0V,0,0,Ct,-St) #Rotate back
	
return(uV,v,w)

end