
#includet(raw""C:\Users\Berlin\.julia\packages\CutAndDisplaceJulia\EW2Bp\test\test_PropagationTest.jl"")

#Our working dir
OuterDir=raw"C:\Users\timmm\Desktop\MeshProp"
cd(OuterDir)
Dir1="WhereTheMeshesLive"
if isdir(Dir1)
	rm(Dir1, recursive=true)
	mkdir(Dir1)
else
	mkdir(Dir1)
end
Dir2="WhereTheResultsLive"
if isdir(Dir2)
	rm(Dir2, recursive=true)
	mkdir(Dir2)
else
	mkdir(Dir2)
end
#Go into mesh dir to start
cd(Dir1)




HSFlag=0; #const
printstyled("Hs off \n",color=:cyan) 

#Elastic constants
G=ShearModulus(2.0e9); 
ν=PoissonsRatio(0.25);#println("Pr is close to 0.5")
(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(G,ν);

#Density
ρrock=2900.;
ρfluid=2600.;
#Gravitational acceleration
g=9.81;

Δρ=((ρrock-ρfluid)*g)

NoTris=300;
errored=[true]
currentNoTris=NoTris
currentKCrit=-79.
CrackVolumeIn=-79.
PropFlag,maxX,minX,maxY,minY,maxZ,minZ=[NaN],[NaN],[NaN],[NaN],[NaN],[NaN],[NaN]


KCrit=1e6; #[5e7 = 50 MPa √m]
for i=10:10:100;

	currentKCrit=KCrit*i

	#Volume
	if Δρ>0
		CrackVolume=((currentKCrit/(Δρ*sqrt(pi)))^(8/3))*((-4*Δρ*(ν-1))/(3*G))
	else
		CrackVolume=((-currentKCrit/(Δρ*sqrt(pi)))^(8/3))*((-4*Δρ*(ν-1))/(3*G))
	end

	for j=1:2
		if j==1
			CrckVolScl=1.62
		else
			CrckVolScl=1.7
		end
		CrackVolumeIn=CrackVolume*CrckVolScl 
		#Running with some catches
		mxlps=4
		x=0
		currentNoTris=NoTris #reset
		while errored[1]==true 
			
			x=x+1
			errored[1]=false
			
			try 
				(PropFlaglp,maxXlp,minXlp,maxYlp,minYlp,maxZlp,minZlp)=testProp(HSFlag,ν,G,Δρ,currentKCrit,CrackVolumeIn,currentNoTris);
				PropFlag[1],maxX[1],minX[1],maxY[1],minY[1],maxZ[1],minZ[1]=PropFlaglp,maxXlp,minXlp,maxYlp,minYlp,maxZlp,minZlp
			catch
				bounds=100 #between -25 and 25
				Pertubation=(round.(rand(1).*bounds))
				currentNoTris=NoTris #reset
				currentNoTris=currentNoTris+Pertubation[1] #try with new no
				errored[1]=true
				printstyled("Errored \n",color=:red) 
			end
			if x==mxlps #spit out the while statement anyway
				errored[1]=false 
				#always errored so set to NaN before write
				PropFlag[1],maxX[1],minX[1],maxY[1],minY[1],maxZ[1],minZ[1]=NaN,NaN,NaN,NaN,NaN,NaN,NaN
			end
			printstyled("$errored[1] \n",color=:red) 
			#if errored[1]==true
			#	error("it errored... why?")
			#end
		end

		#Jump out
		cd(OuterDir)
		println(pwd())
		sleep(2) #julia is too fast
		#remove all meshes
		rm(Dir1, recursive=true)
		#recreate dir
		println("Sleeping 10")
		sleep(10) #stops issues with deleting meshes 
		mkdir(Dir1)

		@info PropFlag[1] maxX[1] minX[1] maxY[1] minY[1] maxZ[1] minZ[1] G ν g Δρ currentNoTris currentKCrit CrackVolumeIn

		cd(Dir2)
		filename="Results-Kcrit-$currentKCrit-Δρ-$Δρ-NoTris-$currentNoTris-GuessAnVolScl-$CrckVolScl"
		OutputDirectory=CutAndDisplaceJulia.LoopResultsWriter(filename,PropFlag[1],maxX[1],minX[1],maxY[1],minY[1],maxZ[1],minZ[1],
										  G,ν,g,Δρ,currentNoTris,currentKCrit,CrackVolumeIn)

		#Jump out
		cd(OuterDir)
		#Back into mesh dir
		cd(Dir1)
		#Reset so the while loop is true again
		errored[1]=true

	end

end

cd(OuterDir)