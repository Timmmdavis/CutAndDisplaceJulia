
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

if Sys.islinux()
	pth='/'
else
	pth='\\'
end	
Dir1FullPath=string(OuterDir,pth,Dir1)
Dir2FullPath=string(OuterDir,pth,Dir2)


Results=readdlm("ResultsTable.txt", ',')

println(Results)

#Go into mesh dir to start
cd(Dir1)

#attempt to stop linux -too many plots error
#GR.inline("png")

HSFlag=0; #const
printstyled("Hs off \n",color=:cyan) 

#Gravitational acceleration
g=9.81;


NoTris=300;
errored=[true]
currentNoTris=NoTris
currentKCrit=-79.
CrackVolumeIn=-79.
PropFlag,maxX,minX,maxY,minY,maxZ,minZ=[NaN],[NaN],[NaN],[NaN],[NaN],[NaN],[NaN]


KCrit=1e6; #[5e7 = 50 MPa √m]
for i=2:size(Results,1)

	#Not bothered if scaling was low - always seems ok
	if Results[i,15]<1.4
		continue
	end
	#Not bothered if it reached the surface
	if Results[i,1]==1
		continue
	end

	#Now set the values
	Δρ=Results[i,11]
	G=Results[i,8]
	ν=Results[i,9]
	currentKCrit=Results[i,13]
	CrckVolScl=Results[i,15]

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

cd(OuterDir)

