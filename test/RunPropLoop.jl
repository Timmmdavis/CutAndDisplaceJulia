
#includet(raw""C:\Users\Berlin\.julia\packages\CutAndDisplaceJulia\EW2Bp\test\test_PropagationTest.jl"")


#Our working dir
OuterDir="/home/tim/Desktop/MeshProp/"
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
if Sys.islinux()
	pth='/'
else
	pth='\\'
end	
Dir1FullPath=string(OuterDir,pth,Dir1)
Dir2FullPath=string(OuterDir,pth,Dir2)

#attempt to stop linux -too many plots error
GR.inline("png")

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
for i=20:20:100;



	for j=1:4 
		if j==1
			CrckVolScl=1.3
		elseif j==2
			CrckVolScl=1.4
		elseif j==3
			CrckVolScl=1.5
		elseif j==4
			CrckVolScl=1.6													
		end					

		for k=1:3

			if k==1
				###########################################
				#Water vs squishy-Granite
				###########################################
				G=ShearModulus(2.0e9); 
				ν=PoissonsRatio(0.25);#println("Pr is close to 0.5")
				(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(G,ν);
				ρrock=2900.;
				ρfluid=1000.;
				Δρ=((ρrock-ρfluid)*g)
				currentKCrit=KCrit*i
			elseif k==2
				###########################################
				#Magma vs Granite
				###########################################
				G=ShearModulus(50.0e9); 
				ν=PoissonsRatio(0.25);#println("Pr is close to 0.5")
				(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(G,ν);
				ρrock=2900.;
				ρfluid=2600.;
				Δρ=((ρrock-ρfluid)*g)
				currentKCrit=KCrit*i				
			elseif k==3
				###########################################
				#Gelatin
				###########################################				
				G=ShearModulus(10e3); 
				ν=PoissonsRatio(0.4999);#println("Pr is close to 0.5")
				(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(G,ν);
				ρrock=1000.;
				ρfluid=1.225;
				Δρ=((ρrock-ρfluid)*g)
				#For gelatin this is really small - 1-10pa
				currentKCrit=i/10 
			end

			#Volume
			if Δρ>0
				CrackVolume=((currentKCrit/(Δρ*sqrt(pi)))^(8/3))*((-4*Δρ*(ν-1))/(3*G))
			else
				CrackVolume=((-currentKCrit/(Δρ*sqrt(pi)))^(8/3))*((-4*Δρ*(ν-1))/(3*G))
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

			#Placing the last png from that loop in the results dir
			filename="Results-Kcrit-$currentKCrit-Δρ-$Δρ-NoTris-$currentNoTris-GuessAnVolScl-$CrckVolScl"
			#Try and move a png IF it exists
			try
				NewPngStr=string(filename,".png")
				ExistingPngStr=CutAndDisplaceJulia.FindingLastPng(Dir1FullPath)
				mv(string(Dir1FullPath,pth,ExistingPngStr),string(Dir2FullPath,pth,NewPngStr); force=true);
				sleep(1) #julia is too fast
			catch
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
			
			OutputDirectory=CutAndDisplaceJulia.LoopResultsWriter(filename,PropFlag[1],maxX[1],minX[1],maxY[1],minY[1],maxZ[1],minZ[1],
											  G,ν,g,Δρ,currentNoTris,currentKCrit,CrackVolumeIn,CrckVolScl)

			#Jump out
			cd(OuterDir)
			#Back into mesh dir
			cd(Dir1)
			#Reset so the while loop is true again
			errored[1]=true

		end

	end

end

cd(OuterDir)