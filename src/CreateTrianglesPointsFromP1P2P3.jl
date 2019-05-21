function CreateTrianglesPointsFromP1P2P3(P1,P2,P3)
#Collapses back for use with external functions in CGAL etc

n=size(P1,1)
#Big long list
Pnts=[P1;P2;P3]

#Reduce and get index's of the new reduced points relative to old index's
Points=unique(Pnts,dims=1) 
PntsUniqueRowNos=UnqiueRowNumberFinder(Pnts,Points)



P1PntsUnique=PntsUniqueRowNos[1:n]
P2PntsUnique=PntsUniqueRowNos[n+1:2*n]
P3PntsUnique=PntsUniqueRowNos[n*2+1:3*n]

#Now place P1 P2 P3 on seperate rows
Points=[1:size(Points,1) Points]
Triangles=[P1PntsUnique P2PntsUnique P3PntsUnique]

return Triangles,Points
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