
#includet(raw""C:\Users\Berlin\.julia\packages\CutAndDisplaceJulia\EW2Bp\test\test_PropagationTest.jl"")


#Our working dir
if Sys.islinux()
	OuterDir="/home/tim/Desktop/MeshProp"
else
	OuterDir=raw"C:\Users\timmm\Desktop\MeshProp"
end

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
#GR.inline("png")

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

NoTris=650;
errored=[true]
currentNoTris=NoTris
currentKCrit=-79.
CrackVolumeIn=-79.
PropFlag,maxX,minX,maxY,minY,maxZ,minZ=[NaN],[NaN],[NaN],[NaN],[NaN],[NaN],[NaN]

list=1:10:101

KCrit=1e6; #[5e7 = 50 MPa √m]
for i=1:length(list)

    #Go from 101:10:1
    sclbackwards=(list[end]+1)-list[i];

	for j=1:10 
		if j==1
			CrckVolScl=2
		elseif j==2
			CrckVolScl=1.9
		elseif j==3
			CrckVolScl=1.8
		elseif j==4
			CrckVolScl=1.7												
		elseif j==5
			CrckVolScl=1.6
		elseif j==6
			CrckVolScl=1.5													
		elseif j==7
			CrckVolScl=1.4
		elseif j==8
			CrckVolScl=1.3												
		elseif j==9
			CrckVolScl=1.2
		elseif j==10
			CrckVolScl=1													
		end				

		for k=1:3

			if k==1
				###########################################
				#Cheb Basin
				###########################################
				ρfluid=469.;
				ρrock=2700.;
				Δρ=((ρfluid-ρrock)*g)
				ν=PoissonsRatio(0.25);#println("Pr is close to 0.5")
				G=ShearModulus(50e9);
				(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(G,ν);
				currentKCrit=KCrit*sclbackwards	
			elseif k==2
				###########################################
				#Pohang
				###########################################
				ρfluid=1000.;
				ρrock=2625.;
				Δρ=((ρfluid-ρrock)*g)
				ν=PoissonsRatio(0.265);#println("Pr is close to 0.5")
				G=ShearModulus(30e9);
				(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(G,ν);
				currentKCrit=KCrit*sclbackwards		
			elseif k==3
				###########################################
				#Gelatin
				###########################################				
				ρfluid=1.225;
				ρrock=1000.;
				Δρ=((ρfluid-ρrock)*g)	
				ν=PoissonsRatio(0.4999);#println("Pr is close to 0.5")
				G=ShearModulus(500);
				(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(G,ν);
				#For gelatin this is really small - 1-100pa
				currentKCrit=sclbackwards
			end

			#Volume (low)
			PKIc2=(currentKCrit^2)*pi
			CrackVolume=real(-(3^(2/3)*PKIc2^(4/3)*(ν - 1))/(16*complex(Δρ)^(5/3)*G));

			CrackVolumeIn=CrackVolume*CrckVolScl 

			#@info CrackVolume Δρ G ν currentKCrit

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
			filename="Results-Kcrit-$currentKCrit-Δρ-$Δρ-NoTris-$currentNoTris-GuessAnVolScl-$CrckVolScl-nu-$ν-mu-$G"
			#Try and move a png IF it exists
			try
				#NewPngStr=string(filename,".png")
				#ExistingPngStr=CutAndDisplaceJulia.FindingLastPng(Dir1FullPath)
				#mv(string(Dir1FullPath,pth,ExistingPngStr),string(Dir2FullPath,pth,NewPngStr); force=true);
				cd(OuterDir)
				mv(Dir1FullPath,string(Dir2FullPath,pth,filename); force=true);
				sleep(1) #julia is too fast #Sleeps 2 mins - should be enough to copy
			catch
			end
			sleep(120) #long enough to copy

			#Jump out
			cd(OuterDir)
			println(pwd())
			sleep(2) #julia is too fast

			#remove all meshes - copied already
			#rm(Dir1, recursive=true)
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