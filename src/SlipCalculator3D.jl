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
	BoundaryConditions=BoundaryConditions.MixedBoundaryConditions
	Stress=BoundaryConditions.Stresses;
	Traction=BoundaryConditions.Tractions;
	#Before calling slip calc check if we will overcome the frictional strength
    (Tn,Tds,Tss)= CalculateNormalAndShearTractions3D( Stress.σxx,Stress.σyy,Stress.σzz,Stress.σxy,Stress.σxz,Stress.σyz,FaceNormalVector );
    #Adding tractions imported into function if these also exist.    
    Tn=Tn+Traction.Tn;
    Tds=Tds+Traction.Tds;
    Tss=Tss+Traction.Tss;
    FricRes=Tn.*µ; #Negative is compression
    if all(-abs.(Tss).>FricRes) && all(-abs.(Tds).>FricRes)
    	println("Cant overcome frictional strength")
        Dn=zeros(size(Tn))
        Dss=zeros(size(Tn))
        Dds=zeros(size(Tn))
        return Dn,Dss,Dds
    end

	#Call default slip calc to get inf matrix and displacements due to BoundaryConditions
	(Dn, Dss, Dds, A, b)=SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)


	#Adding some scaling parameters (Improves Friction Solver performance). 
	#We scale by the average triangle size and the shear mod. 
	(Area,HalfPerimeter ) = AreaOfTriangle3D( P1[:,1],P1[:,2],P1[:,3],P2[:,1],P2[:,2],P2[:,3],P3[:,1],P3[:,2],P3[:,3] );
	Scl=(Statistics.mean(HalfPerimeter)/G); 
	A=A*Scl; 

	#Invert Inf Mat and put inside structs
	InvertedInfMatA=InfMat(inv(A)); 
	b=BoundaryConditionsVec(b);

	#Pass Inf Mat 2 fric func where friction solver is run
	(FricMatPrepped,FricVectorWithoutDisp,L1,L2,L3,L4,L5,Scl,D,
		Vects,Arrys,Flts,Ints,Mats,Bls,IntArrys)=SlipCalculator3D(Scl,n,InvertedInfMatA,b,µ,Sf)

	(Dn,Dss,Dds)=SlipCalculator3D(FricMatPrepped,FricVectorWithoutDisp,L1,L2,L3,L4,L5,Scl,D,
								  Vects,Arrys,Flts,Ints,Mats,Bls,IntArrys)





end


function SlipCalculator3D(Scl,n,InvertedInfMatA::InfMat,b::BoundaryConditionsVec,µ,Sf)
	

	#Each col in [C] represents how much each 
	#element must displace to cause a traction
	#of one unit at element i. Inverted inf mat A

	C=InvertedInfMatA.A #already inverted
	D=C*b.b;

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

	## Construct column vector [b].f
	#b =[D[L1]-CDnTss*Sf-CDnTds*Sf; 
	#	 D[L2]-CDssTss*Sf-CDssTds*Sf;
	#	 D[L3]-CDdsTss*Sf-CDdsTds*Sf;
	#	 2*Sf; 
	#	 2*Sf];

	 # Construct column vector [b].f
	 b =[-CDnTss*Sf-CDnTds*Sf; 
	 	 -CDssTss*Sf-CDssTds*Sf;
	 	 -CDdsTss*Sf-CDdsTds*Sf;
	 	 2*Sf; 
	 	 2*Sf];

	FricMatPrepped=InfMat(a)
	FricVectorWithoutDisp=BoundaryConditionsVec(b)
	#Get arrays ready for fischer_newton
	(Vects,Arrys,Flts,Ints,Mats,Bls,IntArrys)=FischerNewton.InitArrays(n*5);


	return(FricMatPrepped,FricVectorWithoutDisp,L1,L2,L3,L4,L5,Scl,D,Vects,Arrys,Flts,Ints,Mats,Bls,IntArrys)

end

function SlipCalculator3D(FricMatPrepped::InfMat,FricVectorWithoutDisp::BoundaryConditionsVec,L1,L2,L3,L4,L5,Scl,D,
						  Vects,Arrys,Flts,Ints,Mats,Bls,IntArrys)

	a=FricMatPrepped.A
	b=FricVectorWithoutDisp.b
	b[L1] .+=D[L1]; 
	b[L2] .+=D[L2]; 
	b[L3] .+=D[L3]; 

	###################################################################                                                    
	# Solve the linear complementarity problem (calling function that does this)
	###################################################################                                                    
	println("StartingLCP")

	#######path
	#x= pathlcp(ab);
	#######path

	#######LCP solve
	@time x = FischerNewton.fischer_newton(a,b,Vects,Arrys,Flts,Ints,Mats,Bls,IntArrys); 
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


#Friction:
function SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions::MixedBoundaryConditionsFluidVolume,FractureFlag)

	#Compute the tractions acting on the crack
	n=length(FaceNormalVector[:,1]);

	#Extract params and repackage
	Volume=BoundaryConditions.Volumes;
	BoundaryConditions=BoundaryConditions.MixedBoundaryConditions

	FixedEls=zeros(size(FractureFlag));
	#Call default slip calc to get inf matrix and displacements due to BoundaryConditions
	(Dn, Dss, Dds, A, b)=SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)

	#Adding some scaling parameters (Improves solver performance). 
	#We scale by the average triangle size and the shear mod. 
	(Area,HalfPerimeter ) = AreaOfTriangle3D( P1[:,1],P1[:,2],P1[:,3],P2[:,1],P2[:,2],P2[:,3],P3[:,1],P3[:,2],P3[:,3] );
	Scl=(mean(HalfPerimeter)/G); 
	A=A.*Scl;
	#Scl=1;println("Scl off!")  

	#Therefore each col in [C] represents how much each 
	#element must displace to cause a traction
	#of one unit at element i.
	Ainv = inv(A);	

	Norm=maximum(abs.(b[1:n,:])); #Maximum Tn acting on els to close these 
	println("Norm based on max closing stress")
	@info Norm

	AinvF=InfMat(Ainv); 
	println("Setting arbitary fric params")
	µ=fill(0.6,n);
	Sf=zeros(n);
	B=BoundaryConditionsVec(b);
	#Prep our fric mat early (no need to reallocate everytime)
	(FricMatPrepped,FricVectorWithoutDisp,L1,L2,L3,L4,L5,Scl,D,
		Vects,Arrys,Flts,Ints,Mats,Bls,IntArrys)=
	SlipCalculator3D(Scl,n,AinvF,B,µ,Sf)


	NumOfFractures=maximum(FractureFlag);
    NumOfFractures=convert(Int64,NumOfFractures)
    Area=vec(Area)
    #println("SHOULD BE any(FractureFlag.>1) !!!!")
	if any(FractureFlag.>1) #More than two cracks, need sim anneal
        
        X0=zeros(NumOfFractures);
        DesiredAverageHeight=zeros(NumOfFractures);
        Norm=zeros(NumOfFractures);
        DnAvg=zeros(n);
        for i=1:NumOfFractures
        	println("We have this $Volume")
        	println("i is $i")
            if Volume[i]<0 #Neg vols need neg pressure input
                X0[i]=-1.0;  #if two cracks X0=[1,1];
            else
                X0[i]=1.0;  #if two cracks X0=[1,1];
            end

		    #Get the Current crack
		    Indx=FractureFlag.==i;

            AreaFrak_n=0.0
		    for j=eachindex(Indx)
		        if Indx[j]==true
		        	AreaFrak_n=AreaFrak_n+Area[j]
		        end
		    end
            #Compting a good start prssure for each crack given the known opening
            DesiredAverageHeight=(Volume[i]/sum(AreaFrak_n)); error("Check scl here - what do we use for one crack below")#./Scl; 
		    for j=eachindex(Indx)
		        if Indx[j]==true
		            DnAvg[j]=DnAvg[j]+DesiredAverageHeight;
		       end
		    end


        end
		AvgDisps=[DnAvg; zeros(n); zeros(n)]
		ApproxTractions=A*AvgDisps;
		TnApprox=ApproxTractions[1:n];
		maxy=0.0;
		for i=1:NumOfFractures
			Indx=FractureFlag.==i;
			for j=1:length(TnApprox)
				if Indx[j]==true
					if maxy<TnApprox[j]
						maxy=TnApprox[j]
					end
				end
			end
			Norm[i]=maxy;
        end
        Norm=Norm.*10;
		println("Norm based on max traction needed for constant opening desired vol")
		@info Norm        

		X0=vec(X0)
        ## Option 1:
        #Objective function to pass to the simulated annealing solver: 
        ReturnVol=0; #Flag that means the obj func returns volumes 
        ObjectiveFunction=x->ComputePressurisedCrackDn(x,FractureFlag,b,Ainv,Scl,Area,Norm,n,Volume,
        ReturnVol,NumOfFractures,FricMatPrepped,FricVectorWithoutDisp,L1,L2,L3,L4,L5,
        D,Vects,Arrys,Flts,Ints,Mats,Bls,IntArrys)
        #Start value (multiplied by the scalar). 
        #Run the annealing. 
        (InternalPressures,fval) = anneal(ObjectiveFunction, X0);
    	OptimalPressure=Tractions(InternalPressures,[],[]);
    else


		#Compting a good start prssure for each crack given the known opening
		DesiredAverageHeight=(Volume/sum(Area))./Scl; 
		AvgDisps=[ones(n).*DesiredAverageHeight; zeros(n); zeros(n)]
		ApproxTractions=A*AvgDisps; #ApproxTractions
		#ApproxTractions=(A*AvgDisps)/Scl; #ApproxTractions
		#println(maximum(ApproxTractions[1:n]));

		Norm=maximum(ApproxTractions[1:n])
		println("Norm based on max traction needed for constant opening desired vol")
		@info Norm
		#Norm=Norm+1e8

		#Test this norm is going to work: if not fix it...
    	HugeVol=ComputePressurisedCrackDn(10,FractureFlag,b,Ainv,Scl,Area,Norm,n,Volume,
        1,NumOfFractures,FricMatPrepped,FricVectorWithoutDisp,L1,L2,L3,L4,L5,
        D,Vects,Arrys,Flts,Ints,Mats,Bls,IntArrys)
        println("Before")
        @info HugeVol Volume[1] Norm
    	while HugeVol<Volume[1]
    		Norm=Norm*1.1;
    		println("Scaling norm up")
    		HugeVol=ComputePressurisedCrackDn(10,FractureFlag,b,Ainv,Scl,Area,Norm,n,Volume,
	        1,NumOfFractures,FricMatPrepped,FricVectorWithoutDisp,L1,L2,L3,L4,L5,
	        D,Vects,Arrys,Flts,Ints,Mats,Bls,IntArrys)
    	end
    	Reduce=0.05
    	while HugeVol>Volume[1]*10
    		Norm=Norm*(1-Reduce);
    		println("Scaling norm down")
    		HugeVol=ComputePressurisedCrackDn(10,FractureFlag,b,Ainv,Scl,Area,Norm,n,Volume,
	        1,NumOfFractures,FricMatPrepped,FricVectorWithoutDisp,L1,L2,L3,L4,L5,
	        D,Vects,Arrys,Flts,Ints,Mats,Bls,IntArrys)
	        if HugeVol<Volume[1]
	        	Norm=Norm/(1-Reduce);	#Put back and escape
	        	Reduce=Reduce/2; #Reduce size
	        	#Recomp before end of loop
	        	HugeVol=ComputePressurisedCrackDn(10,FractureFlag,b,Ainv,Scl,Area,Norm,n,Volume,
		        1,NumOfFractures,FricMatPrepped,FricVectorWithoutDisp,L1,L2,L3,L4,L5,
		        D,Vects,Arrys,Flts,Ints,Mats,Bls,IntArrys)
	        end
    	end
    	println("After")
    	@info HugeVol Volume[1] Norm

        # Option 2:
        #Objective function to pass to the simple solver: 
        println("rturn vol is 1 for walk and interp")
        ReturnVol=0; #Flag that means the obj func returns volumes 
        ObjectiveFunction=x->ComputePressurisedCrackDn(x,FractureFlag,b,Ainv,Scl,Area,Norm,n,Volume,
        ReturnVol,NumOfFractures,FricMatPrepped,FricVectorWithoutDisp,L1,L2,L3,L4,L5,
        D,Vects,Arrys,Flts,Ints,Mats,Bls,IntArrys)
        #Assuming the Obj func returns volume (for given pressure) ^ just
        #change FIRST output in 'ComputePressurisedCrackDn.m' func above. 
       
        println("Using Optim")
        (res) =Optim.optimize(ObjectiveFunction, 0, 10,method=Brent(),abs_tol=0.001) #
        
        #@info maximum(Ainv[:]) minimum(Ainv[:])

        #Catches error if the max is too small - we increase this and run again
        if abs(abs(Volume[1])-abs(Optim.minimum(res)))>abs(Volume[1])

        	printstyled("Max pressure was too small - also switching to GoldenSection not Brent algo \n",color=:red)
        	(Vects,Arrys,Flts,Ints,Mats,Bls,IntArrys)=FischerNewton.InitArrays(length(L1)*5);println("testing reinit")

        	##Show that 50 is enough:
        	HugeVol=ComputePressurisedCrackDn(50,FractureFlag,b,Ainv,Scl,Area,Norm,n,Volume,
        	        1,NumOfFractures,FricMatPrepped,FricVectorWithoutDisp,L1,L2,L3,L4,L5,
        	        D,Vects,Arrys,Flts,Ints,Mats,Bls,IntArrys)
        	println("Volume with x of 50: $HugeVol")

        	(res) =Optim.optimize(ObjectiveFunction, [9,], Newton(),Optim.Options(iterations = 10,
        																		  store_trace = true,
        																		  show_trace = true,
        																		  f_calls_limit=100)) #
        	#(res) =Optim.optimize(ObjectiveFunction, 0, 50,method=GoldenSection()) #

        end
      	println(summary(res))
      	println(Optim.minimizer(res))
      	println(Optim.minimum(res))
      	println(Optim.iterations(res))
		println(Optim.iteration_limit_reached(res))
      	
       	InternalPressures=Optim.minimizer(res)
       	#(InternalPressures) = WalkAndInterp(ObjectiveFunction, 1e-9, 10, 50,Volume);
        #Put in struct
        OptimalPressure=Tractions(InternalPressures,[],[]);
    end

    #Now compute the min result, more cracks require more output args
    println(OptimalPressure)
    (Dn,Dss,Dds)=ComputePressurisedCrackDn(OptimalPressure,FractureFlag,b,Ainv,Scl,Area,Norm,n,Volume,
        ReturnVol,NumOfFractures,FricMatPrepped,FricVectorWithoutDisp,L1,L2,L3,L4,L5,
        D,Vects,Arrys,Flts,Ints,Mats,Bls,IntArrys)
    
    for i=1:NumOfFractures #For each crack
        Indx=findall(FractureFlag.==i);
        CurrentVol=sum(Dn[Indx].*Area[Indx]);
        #Display results
        println("Volume we want:")
        println(Volume[i])
        println("Volume we have:")
        println(CurrentVol)
         #Catch issues
         if abs(CurrentVol)>abs.(Volume[i]*2)
             println("didnt converge well")
             Dn= Dn*NaN; #nan all values so not used.  
             Dss= Dss*NaN; #nan all values so not used.  
             Dds= Dds*NaN; #nan all values so not used.  
         end     
    end
  
	return Dn,Dss,Dds

end

