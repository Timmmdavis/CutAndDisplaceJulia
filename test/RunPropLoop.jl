
#includet(raw""C:\Users\Berlin\.julia\packages\CutAndDisplaceJulia\EW2Bp\test\test_PropagationTest.jl"")

#Our working dir
OuterDir=raw"C:\Users\Berlin\Desktop\MeshProp"
cd(OuterDir)
Dir1="WhereTheMeshesLive"
if isdir(Dir1)
else
	mkdir(Dir1)
end
Dir2="WhereTheResultsLive"
if isdir(Dir2)
else
mkdir(Dir2)
end
cd(Dir1)




HSFlag=0; #const
printstyled("Hs off \n",color=:cyan) 

#Elastic constants
G=ShearModulus(2.0e9); 
ν=PoissonsRatio(0.25);#println("Pr is close to 0.5")
(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(G,ν);

#Density
ρrock=2900;
ρfluid=2600;
#Gravitational acceleration
g=9.81;

Δρ=((ρrock-ρfluid)*g)


NoTris=300;
errored=[true]

KCrit=1e7; #[5e7 = 50 MPa √m]
for i=1:10:100

	if i>1
		cd(Dir1)
	end
	currentKCrit=KCrit*i

	#Volume
	if Δρ>0
		CrackVolume=((currentKCrit/(Δρ*sqrt(pi)))^(8/3))*((-4*Δρ*(ν-1))/(3*G))
	else
		CrackVolume=((-currentKCrit/(Δρ*sqrt(pi)))^(8/3))*((-4*Δρ*(ν-1))/(3*G))
	end

	for j=1:2
		if j==1
			CrckVolScl=0.5
		else
			CrckVolScl=2
		end
		CrackVolumeIn=CrackVolume*CrckVolScl 
		#Running with some catches
		mxlps=4
		x=0
		while errored[1]==true 
			
			x=x+1
			errored[1]=false
			currentNoTris=NoTris #reset
			try 
				(PropFlag,maxX,minX,maxY,minY,maxZ,minZ)=testProp(HSFlag,ν,G,Δρ,currentKCrit,CrackVolumeIn,currentNoTris);
			catch
				bounds=50 #between -25 and 25
				Pertubation=(round.(rand(1).*bounds)).-bounds
				currentNoTris=currentNoTris+Pertubation[1] #try with new no
				errored[1]=true
				printstyled("Errored \n",color=:red) 
			end
			if x==mxlps #spit out the while statement anyway
				errored[1]=false 
				#always errored so set to NaN before write
				PropFlag,maxX,minX,maxY,minY,maxZ,minZ=NaN,NaN,NaN,NaN,NaN,NaN,NaN
			end
			printstyled("$errored[1] \n",color=:red) 
			if errored[1]==true
				error("it errored... why?")
			end
		end

		#Jump out
		cd(OuterDir)
		#remove all meshes
		rm(Dir1, recursive=true)
		#recreate dir
		mkdir(Dir1)

		cd(Dir2)
		filename="Results-Kcrit-$currentKCrit-Δρ-$Δρ-NoTris-$currentNoTris-GuessAnVolScl-$CrckVolScl"
		OutputDirectory=LoopResultsWriter(filename,PropFlag,maxX,minX,maxY,minY,maxZ,minZ,
										  G.G,ν.ν,g,Δρ,currentNoTris,currentKCrit)

		#Jump out
		cd(OuterDir)

	end

end