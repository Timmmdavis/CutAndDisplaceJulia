function AngSetupStrainFSC(X,Y,Z,bX,bY,bZ,PA,PB,mu,lambda,nu)
# AngSetupFSC_S calculates the Free Surface Correction to strains and 
# stresses associated with angular dislocation pair on each TD side.

# Calculate TD side vector and the angle of the angular dislocation pair
(SideVec,eZ,beta)=CalcSideVec(PA,PB)

if abs(beta)<eps() || abs(pi-beta)<eps()
    Exx = zeros(length(X),1);
	Eyy = zeros(length(X),1);
	Ezz = zeros(length(X),1);
	Exy = zeros(length(X),1);
	Exz = zeros(length(X),1);
	Eyz = zeros(length(X),1);
else
    (b1,b2,b3,I,y1A,y2A,y3A,y1B,y2B,y3B,ey1,ey2,ey3)=CalcSlipVectorDiscCoords(SideVec,eZ,X,Y,Z,PA,bX,bY,bZ,beta)
	
	#AFlip = [ey1[1] ey2[1] ey3[1] ey1[2] ey2[2] ey3[2] ey1[3] ey2[3] ey3[3]]; # Transformation matrix
	#AFlip2 = [ey1 ey2 ey3]; # Transformation matrix
	AFlip = [ey1[1] ey1[2] ey1[3] ey2[1] ey2[2]  ey2[3]  ey3[1] ey3[2]  ey3[3]]; # Transformation matrix
	#@info AFlip AFlip2 AFlip3
	
    Iflp=.!I; #Invert the bool
	indx=findall(I);
	indxf=findall(Iflp);
	
    # For singularities at surface
    v11A = zeros(length(X),1);
    v22A = zeros(length(X),1);
    v33A = zeros(length(X),1);
    v12A = zeros(length(X),1);
    v13A = zeros(length(X),1);
    v23A = zeros(length(X),1);
    
    v11B = zeros(length(X),1);
    v22B = zeros(length(X),1);
    v33B = zeros(length(X),1);
    v12B = zeros(length(X),1);
    v13B = zeros(length(X),1);
    v23B = zeros(length(X),1);
    
	#Init some vars
	b=pi-beta;
	sinB = sin(b);
	cosB = cos(b);
	cotB = cot(b);
	
    # Configuration I
	for i=1:length(indx)
		(v11A[indx[i]],v22A[indx[i]],v33A[indx[i]],v12A[indx[i]],v13A[indx[i]],v23A[indx[i]]) = AngDisStrainFSC(-y1A[indx[i]],-y2A[indx[i]],y3A[indx[i]],cosB,sinB,cotB,-b1,-b2,b3,nu,-PA[3]);
		v13A[indx[i]] = -v13A[indx[i]];
		v23A[indx[i]] = -v23A[indx[i]];
    
		(v11B[indx[i]],v22B[indx[i]],v33B[indx[i]],v12B[indx[i]],v13B[indx[i]],v23B[indx[i]]) = AngDisStrainFSC(-y1B[indx[i]],-y2B[indx[i]],y3B[indx[i]],cosB,sinB,cotB,-b1,-b2,b3,nu,-PB[3]);
		v13B[indx[i]] = -v13B[indx[i]];
		v23B[indx[i]] = -v23B[indx[i]];
	end
    
	b=beta;
	sinB = sin(b);
	cosB = cos(b);
	cotB = cot(b);
    # Configuration II
	for i=1:length(indxf)
		(v11A[indxf[i]],v22A[indxf[i]],v33A[indxf[i]],v12A[indxf[i]],v13A[indxf[i]],v23A[indxf[i]]) = AngDisStrainFSC(y1A[indxf[i]],y2A[indxf[i]],y3A[indxf[i]],cosB,sinB,cotB,b1,b2,b3,nu,-PA[3]);
		
		(v11B[indxf[i]],v22B[indxf[i]],v33B[indxf[i]],v12B[indxf[i]],v13B[indxf[i]],v23B[indxf[i]]) = AngDisStrainFSC(y1B[indxf[i]],y2B[indxf[i]],y3B[indxf[i]],cosB,sinB,cotB,b1,b2,b3,nu,-PB[3]);
	end



    # Calculate total Free Surface Correction to strains in ADCS
    v11 = v11B-v11A;
    v22 = v22B-v22A;
    v33 = v33B-v33A;
    v12 = v12B-v12A;
    v13 = v13B-v13A;
    v23 = v23B-v23A;
	

    # Transform total Free Surface Correction to strains from ADCS to EFCS
    (Exx,Eyy,Ezz,Exy,Exz,Eyz) = TensorTransformation3D(v11,v22,v33,v12,v13,v23,AFlip);
end

return(Exx,Eyy,Ezz,Exy,Exz,Eyz);
	
end