function SlipCalculator3D(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,BoundaryConditions,FixedEls)

	#Compute the tractions acting on the crack
	n=length(FaceNormalVector[:,1]);
	(Tn,Tds,Tss)=SetupTractionVector(BoundaryConditions,FaceNormalVector,n,G,λ);

	(TractionInfMats)=CalculateInfluenceMatrices3D(FaceNormalVector,MidPoint,P1,P2,P3,ν,G,λ,FixedEls,HSFlag,n)

	println("Put me inside SetupTractionVector")
	FixedFlag=FixedEls.==1
	NotFixedFlag=FixedEls.!=1
	if any(FixedFlag)
		#Computing for fixed els
		for j=1:length(FixedFlag)
			if FixedFlag[j]==true
				Tn[j]=0.0;
				Tds[j]=0.0;
				Tss[j]=0.0;
			end
		end
	end

	#Concatenate influence matrix ready for equation 
	Atn  = [TractionInfMats.DnTn  TractionInfMats.DssTn  TractionInfMats.DdsTn ];     
	Atss = [TractionInfMats.DnTss TractionInfMats.DssTss TractionInfMats.DdsTss];     
	Atds = [TractionInfMats.DnTds TractionInfMats.DssTds TractionInfMats.DdsTds];     
	A= -[Atn;Atss;Atds];  

	#Prep traction vector
	b=[Tn; Tss; Tds];

	#Do linear equation
	D=A\b;


	n=sum(NotFixedFlag)
	#Extract arrays
	Dn=D[1:n];
	Dss=D[n+1:2*n];
	Dds=D[n*2+1:3*n];

	return Dn, Dss, Dds

end

