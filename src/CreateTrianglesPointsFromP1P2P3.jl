function CreateTrianglesPointsFromP1P2P3(P1,P2,P3)
#Collapses back for use with external functions in CGAL etc




#Concat 2 big mat
AllPnts=CreateSortedPointsP1P2P3(P1,P2,P3);
println(size(AllPnts))
println(size(P1))
println(P1[end,:])
println(P2[end,:])
println(P3[end,:])
println(AllPnts[end,:])
error("This should be sorted!")

Slimmed=unique(AllPnts,dims=1) 

#Extract as these are now sorted:
P1=copy(AllPnts[:,1:3])
P2=copy(AllPnts[:,4:6])
P3=copy(AllPnts[:,7:9])

#Reduce and get index's of the new reduced points relative to old index's
P1Unique=unique(P1,dims=1) 
P1RowNos=UnqiueRowNumberFinder(P1,P1Unique)

#Reduce and get index's of the new reduced points relative to old index's
P2Unique=unique(P2,dims=1) 
P2RowNos=UnqiueRowNumberFinder(P2,P2Unique)

#Reduce and get index's of the new reduced points relative to old index's
P3Unique=unique(P3,dims=1) 
P3RowNos=UnqiueRowNumberFinder(P3,P1Unique)


#Now place P1 P2 P3 on seperate rows
n=length(P1Unique)
Points=[1:n P1Unique P2Unique P3Unique]
Triangles=[P1RowNos P2RowNos P3RowNos]

#=
Points=zeros(Int(length(P1)),3)
Points[1:3:end,:]=copy(P1);
Points[2:3:end,:]=copy(P2);
Points[3:3:end,:]=copy(P3);
n=length(Points)
Points=[1:n/3 Points]
Triangles=fill(0,Int(n/9),3)
Triangles[:,1]=1:3:n/3;
Triangles[:,2]=2:3:n/3;
Triangles[:,3]=3:3:n/3;
=#

return Points,Triangles
end 



function UnqiueRowNumberFinder(OriginalData,UnqiueData)

n=size(OriginalData,1)
n2=size(OriginalData,2)
RowNos=zeros(n)
IndividualRows=fill(false,n)

for i=1:n2
	RowNos=indexin(OriginalData[:,i],UnqiueData[:,i])
	for j=1:length(RowNos)
		if RowNos[j]!=j
			IndividualRows[j]=true
		end
	end
end
for j=1:length(RowNos)
	if IndividualRows[j]==true
		RowNos[j]=j
	end
end
RowNos=Int.(RowNos)
end