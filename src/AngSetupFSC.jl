function AngSetupFSC(X,Y,Z,bX,bY,bZ,PA,PB,nu)
# AngSetupFSC calculates the Free Surface Correction to displacements 
# associated with angular dislocation pair on each TD side.

# Calculate TD side vector and the angle of the angular dislocation pair
SideVec = PB-PA;
eZ = [0;0;1];

G=-SideVec'*eZ/norm(SideVec);
beta = acos(G[1]);

if abs(beta)<eps() || abs(pi-beta)<eps()
    ue = zeros(length(X),1);
    un = zeros(length(X),1);
    uv = zeros(length(X),1);
else
    ey1 = [SideVec[1:2];0];
    ey1 = ey1/norm(ey1);
    ey3 = -eZ;
    ey2 = cross(ey3,ey1);
	
    # Transform coordinates from EFCS to the first ADCS
	(y1A,y2A,y3A)=RotateObject3DNewCoords(X,Y,Z,PA[1],PA[2],PA[3],ey1,ey2,ey3)
    
	# Transform coordinates from EFCS to the second ADCS
	(y1AB,y2AB,y3AB)=RotateObject3DNewCoords(SideVec[1],SideVec[2],SideVec[3],0,0,0,ey1,ey2,ey3)
    y1B = y1A.-y1AB;
    y2B = y2A.-y2AB;
    y3B = y3A.-y3AB;
    
    # Transform slip vector components from EFCS to ADCS
	(b1,b2,b3)=RotateObject3DNewCoords(bX,bY,bZ,0,0,0,ey1,ey2,ey3)
    
    # Determine the best arteact-free configuration for the calculation
    # points near the free furface
    I = (beta.*y1A).>=0;
    
	#InitOutputs
	v1A= Array{Float64}(undef, length(y1A),1);
	v2A= Array{Float64}(undef, length(y1A),1);
	v3A= Array{Float64}(undef, length(y1A),1);
	v1B= Array{Float64}(undef, length(y1A),1);
	v2B= Array{Float64}(undef, length(y1A),1);
	v3B= Array{Float64}(undef, length(y1A),1);
	
	Iflp=.!I; #Invert the bool

	indx=findall(I);
	indxf=findall(Iflp);
	b=-pi+beta;
	sinB = sin(b);
	cosB = cos(b);
	cotB = cot(b);
	cotB2=cot(b/2);
	# Configuration I
	for i=1:length(indx)
		(v1A[indx[i]],v2A[indx[i]],v3A[indx[i]]) = AngDisDispFSC(y1A[indx[i]],y2A[indx[i]],y3A[indx[i]],cosB,sinB,cotB,cotB2,b1,b2,b3,nu,-PA[3]);
		(v1B[indx[i]],v2B[indx[i]],v3B[indx[i]]) = AngDisDispFSC(y1B[indx[i]],y2B[indx[i]],y3B[indx[i]],cosB,sinB,cotB,cotB2,b1,b2,b3,nu,-PB[3]);
	end
	b=beta;
	sinB = sin(b);
	cosB = cos(b);
	cotB = cot(b);
	cotB2=cot(b/2);
	# Configuration II
	for i=1:length(indxf)
		(v1A[indxf[i]],v2A[indxf[i]],v3A[indxf[i]]) = AngDisDispFSC(y1A[indxf[i]],y2A[indxf[i]],y3A[indxf[i]],cosB,sinB,cotB,cotB2,b1,b2,b3,nu,-PA[3]);
		(v1B[indxf[i]],v2B[indxf[i]],v3B[indxf[i]]) = AngDisDispFSC(y1B[indxf[i]],y2B[indxf[i]],y3B[indxf[i]],cosB,sinB,cotB,cotB2,b1,b2,b3,nu,-PB[3]);
    end
	
    # Calculate total Free Surface Correction to displacements in ADCS
    v1 = v1B-v1A;
    v2 = v2B-v2A;
    v3 = v3B-v3A;
    
	#Inverse
	Vx=[ey1[1],ey2[1],ey3[1]];
	Vy=[ey1[2],ey2[2],ey3[2]];
	Vz=[ey1[3],ey2[3],ey3[3]];
	
    # Transform total Free Surface Correction to displacements from ADCS 
    # to EFCS
	(ue,un,uv)=RotateObject3DNewCoords(v1,v2,v3,0,0,0,Vx,Vy,Vz)

end	
return(ue,un,uv)
end