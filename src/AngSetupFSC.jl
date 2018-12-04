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
	
    # Configuration I
    (v1A[I],v2A[I],v3A[I]) = AngDisDispFSC(y1A[I],y2A[I],y3A[I],-pi+beta,b1,b2,b3,nu,-PA[3]);
    (v1B[I],v2B[I],v3B[I]) = AngDisDispFSC(y1B[I],y2B[I],y3B[I],-pi+beta,b1,b2,b3,nu,-PB[3]);
	
	Iflp=.!I; #Invert the bool
	
    # Configuration II
    (v1A[Iflp],v2A[Iflp],v3A[Iflp]) = AngDisDispFSC(y1A[Iflp],y2A[Iflp],y3A[Iflp],beta,b1,b2,b3,nu,-PA[3]);
    (v1B[Iflp],v2B[Iflp],v3B[Iflp]) = AngDisDispFSC(y1B[Iflp],y2B[Iflp],y3B[Iflp],beta,b1,b2,b3,nu,-PB[3]);
    
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