function AngSetupDispFSC(X,Y,Z,bX,bY,bZ,PA,PB,nu,Vnorm,Vstrike,Vdip)
# AngSetupFSC calculates the Free Surface Correction to displacements 
# associated with angular dislocation pair on each TD side.

# Calculate TD side vector and the angle of the angular dislocation pair
(SideVec,eZ,beta)=CalcSideVec(PA,PB)

if abs(beta)<eps() || abs(pi-beta)<eps()
    ue = zeros(length(X),1);
    un = zeros(length(X),1);
    uv = zeros(length(X),1);
else
    (b1,b2,b3,I,y1A,y2A,y3A,y1B,y2B,y3B,ey1,ey2,ey3)=CalcSlipVectorDiscCoords(SideVec,eZ,X,Y,Z,PA,bX,bY,bZ,beta)
    
	println("Make sure there is no allocation in these")
	## Transform slip vector components from TDCS into EFCS
	(Dn1__,Dss0n_,Dds0n_) = RotateObject3DNewCoords(1.,0.,0.,0,0,0,Vnorm,Vstrike,Vdip);
	(Dn0ss,Dss1__,Dds1ss) = RotateObject3DNewCoords(0.,1.,0.,0,0,0,Vnorm,Vstrike,Vdip);
	(Dn0ds,Dss0ds,Dds1__) = RotateObject3DNewCoords(0.,0.,1.,0,0,0,Vnorm,Vstrike,Vdip);
	# Transform slip vector components from EFCS to ADCS
	(Dn1__,Dss0n_,Dds0n_)=RotateObject3DNewCoords(Dn1__,Dss0n_,Dds0n_,0,0,0,ey1,ey2,ey3)
	(Dn0ss,Dss1__,Dds1ss)=RotateObject3DNewCoords(Dn0ss,Dss1__,Dds1ss,0,0,0,ey1,ey2,ey3)
	(Dn0ds,Dss0ds,Dds1__)=RotateObject3DNewCoords(Dn0ds,Dss0ds,Dds1__,0,0,0,ey1,ey2,ey3)
		
	println("Avoid this allocation (coming from above, be wary changing this from below...)")
	#InitOutputs
	Ux  = Array{Float64}(undef, length(X),1);
	Vx  = Array{Float64}(undef, length(X),1);
	Wx  = Array{Float64}(undef, length(X),1);
	Uy  = Array{Float64}(undef, length(X),1);
	Vy  = Array{Float64}(undef, length(X),1);
	Wy  = Array{Float64}(undef, length(X),1);
	Uz  = Array{Float64}(undef, length(X),1);
	Vz  = Array{Float64}(undef, length(X),1);
	Wz  = Array{Float64}(undef, length(X),1);

	
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
	
		(uxA,uyA,uzA,vxA,vyA,vzA,wxA,wyA,wzA) = AngDisDispFSC(y1A[indx[i]],y2A[indx[i]],y3A[indx[i]],cosB,sinB,cotB,cotB2,b1,b2,b3,nu,-PA[3]);
		
		
		#Add mixed components together (different coords)
		Ux[indx[i]] = -((Dn1__/4/pi/(1-nu)*uxA)+(Dss0n_/4/pi/(1-nu)*uyA)+(Dds0n_/4/pi/(1-nu)*uzA))
		Uy[indx[i]] = -((Dn0ss/4/pi/(1-nu)*uxA)+(Dss1__/4/pi/(1-nu)*uyA)+(Dds1ss/4/pi/(1-nu)*uzA))
		Uz[indx[i]] = -((Dn0ds/4/pi/(1-nu)*uxA)+(Dss0ds/4/pi/(1-nu)*uyA)+(Dds1__/4/pi/(1-nu)*uzA))
		Vx[indx[i]] = -((Dn1__/4/pi/(1-nu)*vxA)+(Dss0n_/4/pi/(1-nu)*vyA)+(Dds0n_/4/pi/(1-nu)*vzA))
		Vy[indx[i]] = -((Dn0ss/4/pi/(1-nu)*vxA)+(Dss1__/4/pi/(1-nu)*vyA)+(Dds1ss/4/pi/(1-nu)*vzA))
		Vz[indx[i]] = -((Dn0ds/4/pi/(1-nu)*vxA)+(Dss0ds/4/pi/(1-nu)*vyA)+(Dds1__/4/pi/(1-nu)*vzA))
		Wx[indx[i]] = -((Dn1__/4/pi/(1-nu)*wxA)+(Dss0n_/4/pi/(1-nu)*wyA)+(Dds0n_/4/pi/(1-nu)*wzA))
		Wy[indx[i]] = -((Dn0ss/4/pi/(1-nu)*wxA)+(Dss1__/4/pi/(1-nu)*wyA)+(Dds1ss/4/pi/(1-nu)*wzA))
		Wz[indx[i]] = -((Dn0ds/4/pi/(1-nu)*wxA)+(Dss0ds/4/pi/(1-nu)*wyA)+(Dds1__/4/pi/(1-nu)*wzA))
		
		
		(uxB,uyB,uzB,vxB,vyB,vzB,wxB,wyB,wzB) = AngDisDispFSC(y1B[indx[i]],y2B[indx[i]],y3B[indx[i]],cosB,sinB,cotB,cotB2,b1,b2,b3,nu,-PB[3]);
		
		#Add mixed components together (different coords)
		Ux[indx[i]] =Ux[indx[i]] + ((Dn1__/4/pi/(1-nu)*uxB)+(Dss0n_/4/pi/(1-nu)*uyB)+(Dds0n_/4/pi/(1-nu)*uzB))
		Uy[indx[i]] =Uy[indx[i]] + ((Dn0ss/4/pi/(1-nu)*uxB)+(Dss1__/4/pi/(1-nu)*uyB)+(Dds1ss/4/pi/(1-nu)*uzB))
		Uz[indx[i]] =Uz[indx[i]] + ((Dn0ds/4/pi/(1-nu)*uxB)+(Dss0ds/4/pi/(1-nu)*uyB)+(Dds1__/4/pi/(1-nu)*uzB))
		Vx[indx[i]] =Vx[indx[i]] + ((Dn1__/4/pi/(1-nu)*vxB)+(Dss0n_/4/pi/(1-nu)*vyB)+(Dds0n_/4/pi/(1-nu)*vzB))
		Vy[indx[i]] =Vy[indx[i]] + ((Dn0ss/4/pi/(1-nu)*vxB)+(Dss1__/4/pi/(1-nu)*vyB)+(Dds1ss/4/pi/(1-nu)*vzB))
		Vz[indx[i]] =Vz[indx[i]] + ((Dn0ds/4/pi/(1-nu)*vxB)+(Dss0ds/4/pi/(1-nu)*vyB)+(Dds1__/4/pi/(1-nu)*vzB))
		Wx[indx[i]] =Wx[indx[i]] + ((Dn1__/4/pi/(1-nu)*wxB)+(Dss0n_/4/pi/(1-nu)*wyB)+(Dds0n_/4/pi/(1-nu)*wzB))
		Wy[indx[i]] =Wy[indx[i]] + ((Dn0ss/4/pi/(1-nu)*wxB)+(Dss1__/4/pi/(1-nu)*wyB)+(Dds1ss/4/pi/(1-nu)*wzB))
		Wz[indx[i]] =Wz[indx[i]] + ((Dn0ds/4/pi/(1-nu)*wxB)+(Dss0ds/4/pi/(1-nu)*wyB)+(Dds1__/4/pi/(1-nu)*wzB))
		
	end
	
	
	b=beta;
	sinB = sin(b);
	cosB = cos(b);
	cotB = cot(b);
	cotB2=cot(b/2);
	# Configuration II
	for i=1:length(indxf)
	
		(uxA,uyA,uzA,vxA,vyA,vzA,wxA,wyA,wzA) = AngDisDispFSC(y1A[indxf[i]],y2A[indxf[i]],y3A[indxf[i]],cosB,sinB,cotB,cotB2,b1,b2,b3,nu,-PA[3]);
		
		#Add mixed components together (different coords)
		Ux[indxf[i]] = -((Dn1__/4/pi/(1-nu)*uxA)+(Dss0n_/4/pi/(1-nu)*uyA)+(Dds0n_/4/pi/(1-nu)*uzA))
		Uy[indxf[i]] = -((Dn0ss/4/pi/(1-nu)*uxA)+(Dss1__/4/pi/(1-nu)*uyA)+(Dds1ss/4/pi/(1-nu)*uzA))
		Uz[indxf[i]] = -((Dn0ds/4/pi/(1-nu)*uxA)+(Dss0ds/4/pi/(1-nu)*uyA)+(Dds1__/4/pi/(1-nu)*uzA))
		Vx[indxf[i]] = -((Dn1__/4/pi/(1-nu)*vxA)+(Dss0n_/4/pi/(1-nu)*vyA)+(Dds0n_/4/pi/(1-nu)*vzA))
		Vy[indxf[i]] = -((Dn0ss/4/pi/(1-nu)*vxA)+(Dss1__/4/pi/(1-nu)*vyA)+(Dds1ss/4/pi/(1-nu)*vzA))
		Vz[indxf[i]] = -((Dn0ds/4/pi/(1-nu)*vxA)+(Dss0ds/4/pi/(1-nu)*vyA)+(Dds1__/4/pi/(1-nu)*vzA))
		Wx[indxf[i]] = -((Dn1__/4/pi/(1-nu)*wxA)+(Dss0n_/4/pi/(1-nu)*wyA)+(Dds0n_/4/pi/(1-nu)*wzA))
		Wy[indxf[i]] = -((Dn0ss/4/pi/(1-nu)*wxA)+(Dss1__/4/pi/(1-nu)*wyA)+(Dds1ss/4/pi/(1-nu)*wzA))
		Wz[indxf[i]] = -((Dn0ds/4/pi/(1-nu)*wxA)+(Dss0ds/4/pi/(1-nu)*wyA)+(Dds1__/4/pi/(1-nu)*wzA))	
		
		(uxB,uyB,uzB,vxB,vyB,vzB,wxB,wyB,wzB) = AngDisDispFSC(y1B[indxf[i]],y2B[indxf[i]],y3B[indxf[i]],cosB,sinB,cotB,cotB2,b1,b2,b3,nu,-PB[3]);
		
		#Add mixed components together (different coords)
		Ux[indxf[i]] =Ux[indxf[i]] + ((Dn1__/4/pi/(1-nu)*uxB)+(Dss0n_/4/pi/(1-nu)*uyB)+(Dds0n_/4/pi/(1-nu)*uzB))
		Uy[indxf[i]] =Uy[indxf[i]] + ((Dn0ss/4/pi/(1-nu)*uxB)+(Dss1__/4/pi/(1-nu)*uyB)+(Dds1ss/4/pi/(1-nu)*uzB))
		Uz[indxf[i]] =Uz[indxf[i]] + ((Dn0ds/4/pi/(1-nu)*uxB)+(Dss0ds/4/pi/(1-nu)*uyB)+(Dds1__/4/pi/(1-nu)*uzB))
		Vx[indxf[i]] =Vx[indxf[i]] + ((Dn1__/4/pi/(1-nu)*vxB)+(Dss0n_/4/pi/(1-nu)*vyB)+(Dds0n_/4/pi/(1-nu)*vzB))
		Vy[indxf[i]] =Vy[indxf[i]] + ((Dn0ss/4/pi/(1-nu)*vxB)+(Dss1__/4/pi/(1-nu)*vyB)+(Dds1ss/4/pi/(1-nu)*vzB))
		Vz[indxf[i]] =Vz[indxf[i]] + ((Dn0ds/4/pi/(1-nu)*vxB)+(Dss0ds/4/pi/(1-nu)*vyB)+(Dds1__/4/pi/(1-nu)*vzB))
		Wx[indxf[i]] =Wx[indxf[i]] + ((Dn1__/4/pi/(1-nu)*wxB)+(Dss0n_/4/pi/(1-nu)*wyB)+(Dds0n_/4/pi/(1-nu)*wzB))
		Wy[indxf[i]] =Wy[indxf[i]] + ((Dn0ss/4/pi/(1-nu)*wxB)+(Dss1__/4/pi/(1-nu)*wyB)+(Dds1ss/4/pi/(1-nu)*wzB))
		Wz[indxf[i]] =Wz[indxf[i]] + ((Dn0ds/4/pi/(1-nu)*wxB)+(Dss0ds/4/pi/(1-nu)*wyB)+(Dds1__/4/pi/(1-nu)*wzB))
		
    end
	

	
	#Inverse rot mat
	VxR=[ey1[1],ey2[1],ey3[1]];
	VyR=[ey1[2],ey2[2],ey3[2]];
	VzR=[ey1[3],ey2[3],ey3[3]];
	
    # Transform total Free Surface Correction to displacements from ADCS 
    # to EFCS
	println("Reduce allocation here")
	(Ux,Vx,Wx)=RotateObject3DNewCoords(Ux,Vx,Wx,0,0,0,VxR,VyR,VzR)
	(Uy,Vy,Wy)=RotateObject3DNewCoords(Uy,Vy,Wy,0,0,0,VxR,VyR,VzR)
	(Uz,Vz,Wz)=RotateObject3DNewCoords(Uz,Vz,Wz,0,0,0,VxR,VyR,VzR)


end	#if statement

return(Ux,Vx,Wx,Uy,Vy,Wy,Uz,Vz,Wz)
end