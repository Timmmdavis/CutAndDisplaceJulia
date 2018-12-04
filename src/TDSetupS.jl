function TDSetupS(x,y,z,alpha,bx,by,bz,nu,TriVertex,SideVec)
# TDSetupS transforms coordinates of the calculation points as well as 
# slip vector components from ADCS into TDCS. It then calculates the 
# strains in ADCS and transforms them into TDCS.

(Ct,St,y1,z1,by1,bz1)=TransformToADCS(y,z,by,bz,SideVec,TriVertex)

#Init arrays
Ang=-pi+alpha;
cosA = cos(Ang);
sinA = sin(Ang);
exxV = Array{Float64}(undef, length(x),1);
eyyV = Array{Float64}(undef, length(x),1);
ezzV = Array{Float64}(undef, length(x),1);
exyV = Array{Float64}(undef, length(x),1);
exzV = Array{Float64}(undef, length(x),1);
eyzV = Array{Float64}(undef, length(x),1);
#Extra defs out of loop to speed it up
E1=(1-nu); #Elastic cons
E2=(2*nu+1);
E3=bx/8/pi/E1;
cosA2=cosA^2;
sinADE1=sinA/8/pi/(1-nu);


#@info "Some variables" x[1] y[1] z[1] cosA sinA bx by1[1] bz1[1] nu E1 cosA2 sinADE1
#poop


# Calculate strains associated with an angular dislocation in ADCS
for i=1:length(x)
	(exx,eyy,ezz,exy,exz,eyz) = AngDisStrain(x[i],y1[i],z1[i],cosA,sinA,bx,by1[1],bz1[1],nu,E1,E2,E3,cosA2,sinADE1)
	exxV[i]=exx;
	eyyV[i]=eyy;
	ezzV[i]=ezz;
	exyV[i]=exy;
	exzV[i]=exz;
	eyzV[i]=eyz;
end	

# Transform strains from ADCS into TDCS
B=[[1 0 0];[0 Ct St];[0 -St Ct]]; # 3x3 Transformation matrix
(exxV,eyyV,ezzV,exyV,exzV,eyzV) = TensorTransformation3D(exxV,eyyV,ezzV,exyV,exzV,eyzV,B);

return(exxV,eyyV,ezzV,exyV,exzV,eyzV)

end