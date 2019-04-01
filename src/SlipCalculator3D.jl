function SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)

	#Compute the tractions acting on the crack
	n=length(FaceNormalVector[:,1]);
	(Tn,Tds,Tss)=SetupTractionVector(BoundaryConditions,FaceNormalVector,n,G,λ);
	(Tn,Tds,Tss)=KnockOutFixedRows(FixedEls,Tn,Tds,Tss)

	#Compute the influence matricies
	(TractionInfMats)=CalculateInfluenceMatrices3D(FaceNormalVector,MidPoint,P1,P2,P3,ν,G,λ,FixedEls,HSFlag,n)
	println("OutOfTD")	

	#Concatenate influence matrix ready for equation 
	Atn  = [TractionInfMats.DnTn  TractionInfMats.DssTn  TractionInfMats.DdsTn ];     
	Atss = [TractionInfMats.DnTss TractionInfMats.DssTss TractionInfMats.DdsTss];     
	Atds = [TractionInfMats.DnTds TractionInfMats.DssTds TractionInfMats.DdsTds];     
	A= -[Atn;Atss;Atds]; 

	#Prep traction vector
	b=[Tn; Tss; Tds];

	#Do linear equation
	D=A\b;
	println("LinearEqDone")	

	n=sum(FixedEls.==0)
	#Extract arrays
	Dn=D[1:n];
	Dss=D[n+1:2*n];
	Dds=D[n*2+1:3*n];

	return Dn, Dss, Dds, A, b

end

#Friction:
function SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions::MixedBoundaryConditionsFriction,FixedEls)

	#Compute the tractions acting on the crack
	n=length(FaceNormalVector[:,1]);

	#Extract params and repackage
	µ=BoundaryConditions.FrictionParameters.µ;
	Sf=BoundaryConditions.FrictionParameters.Sf;
	Stress=BoundaryConditions.Stresses;
	Traction=BoundaryConditions.Tractions;
	BoundaryConditions=MixedBoundaryConditions(Stress,Traction);
	#Call default slip calc to get inf matrix and displacements due to BoundaryConditions
	(Dn, Dss, Dds, A, b)=SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)

	#Adding some scaling parameters (Improves Friction Solver performance). 
	#We scale by the average triangle size and the shear mod. 
	(Area,HalfPerimeter ) = AreaOfTriangle3D( P1[:,1],P1[:,2],P1[:,3],P2[:,1],P2[:,2],P2[:,3],P3[:,1],P3[:,2],P3[:,3] );
	Scl=(mean(HalfPerimeter)/G); 
	A=A.*Scl;  

	#Therefore each col in [C] represents how much each 
	#element must displace to cause a traction
	#of one unit at element i.
	C = inv(A);	

	D=C*b;

	# Construct a and b for the equation [y]=[a][x]+[b] where
	#  [y]  = |+Dn|     and [x] = |-tn|   
	#         |+Dss|              |-tss|   
	#         |+Dds|              |-tds|    
	#         |+Tss|              |-Dss|   
	#         |+Tds|              |-Dds|   
	# Where plus and minus follow the convention that positive (+) represents: 
	# Dn = Opening
	# Dss= Left lat movement
	# Dds= Reverse slip when normals face up
	# Tn = Tension
	# Tss= Positive when counter clockwise from normal
	# Tds= Positive when facing in positive z dir on tris normal side. 

	# Creating some lengths we can use for extraction of sub matricies/vectors
	# n being the number of elements in one edge a submatrix/vector
	L1=1:n;
	L2=n+1:2*n;
	L3=2*n+1:3*n;
	L4=3*n+1:4*n;
	L5=4*n+1:5*n;

	# Now extract variables from inverted matrix [C] 
	#Submatrix in 'column' 1 disps cause traction Tn
	CDnTn   =C[L1,L1];
	CDssTn  =C[L2,L1];	
	CDdsTn  =C[L3,L1];
	#Submatrix in 'column' 2 disps cause traction Tss
	CDnTss  =C[L1,L2];	
	CDssTss =C[L2,L2];	
	CDdsTss =C[L3,L2];
	#Submatrix in 'column' 3 disps cause traction Tds
	CDnTds  =C[L1,L3];  		
	CDssTds =C[L2,L3];	
	CDdsTds =C[L3,L3];	
 
	#Making sure input frictions are vectors
	µ=zeros(n).+µ; 
	Sf=zeros(n).+Sf; 

	# Form n by n array with coefficients of friction on the diagonal.
	dµ = diagm(0 => µ);
	# Allocate n by n identity and zero matrices.#I uses LinearAlgebra package
	ID = Matrix{Float64}(I,n,n) #ID = eye(n);
	ZE = zeros(n,n);

	# Construct matrix [a]. Modified form of Equation 28 - Kaven 2012
	a =	[(CDnTn- CDnTss*dµ-CDnTds*dµ)   CDnTss  CDnTds  ZE ZE;
		 (CDssTn-CDssTss*dµ-CDssTds*dµ) CDssTss CDssTds ID ZE;
		 (CDdsTn-CDdsTss*dµ-CDdsTds*dµ) CDdsTss CDdsTds ZE ID;
		 (2*dµ)                          -ID      ZE     ZE ZE;
		 (2*dµ)                           ZE     -ID     ZE ZE];

	# Construct column vector [b].
	b =	[D[L1]-CDnTss*Sf-CDnTds*Sf; 
		 D[L2]-CDssTss*Sf-CDssTds*Sf;
		 D[L3]-CDdsTss*Sf-CDdsTds*Sf;
		 2*Sf; 
		 2*Sf];

	###################################################################                                                    
	# Solve the linear complementarity problem (calling function that does this)
	###################################################################                                                    
	println("StartingLCP")

	#######path
	#x= pathlcp(ab);
	#######path

	#######LCP solve
	@time x = FischerNewton.fischer_newton(a,b); 
	#######LCP solve

	println("LcpSpeed[s])");
	###################################################################

	#Calculating [y] comp eq: [y]=[a][x]+[b]
	y = a*x+b;
                                                    
	#Extracting sub parts of vectors: 
	#[x]
	x1=x[L1,1];
	x2=x[L2,1];
	x3=x[L3,1];
	x4=x[L4,1];
	x5=x[L5,1];
	#[y]
	y1=y[L1,1];
	y2=y[L2,1];
	y3=y[L3,1];
	y4=y[L4,1];
	y5=y[L5,1];

	#Extracting slip                                      
	Dn  = y1;                            
	Dss = y2-x4;     
	Dds = y3-x5;

	#Extracting the resultant traction at each elements midpoint
	Tn=-x1;
	Tss=y4-x2;
	Tds=y5-x3;
	#This the remote stress left over once stress induced by the slip
	#from the boundary has been removed

	#Scaling back the data. 
	Dn=Dn.*Scl;
	Dss=Dss.*Scl;
	Dds=Dds.*Scl;

	return Dn,Dss,Dds

end

