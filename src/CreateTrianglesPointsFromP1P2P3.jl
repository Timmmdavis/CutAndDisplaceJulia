function CreateTrianglesPointsFromP1P2P3(P1,P2,P3)
#Collapses back for use with external functions in CGAL etc

n=size(P1,1)
#Big long list
Pnts=[P1;P2;P3]

#Reduce and get index's of the new reduced points relative to old index's
Points=unique(Pnts,dims=1) 
PntsUniqueRowNos=UniqueRowNumberFinder(Pnts,Points)
P1PntsUnique=PntsUniqueRowNos[1:n]
P2PntsUnique=PntsUniqueRowNos[n+1:2*n]
P3PntsUnique=PntsUniqueRowNos[n*2+1:3*n]

#Now place P1 P2 P3 on seperate rows
Points=[1:size(Points,1) Points]
Triangles=[P1PntsUnique P2PntsUnique P3PntsUnique]

return Triangles,Points
end 



function UniqueRowNumberFinder(OriginalData,UniqueData)

n=size(OriginalData,1)
m=size(UniqueData,1)

RowNos=fill(0,n)
#Lp to check that each col is unique - if so we then use the index of the 
#new unqiue row
counter=0
flag=true
p=1:3
for i=1:n
	for j=1:m
		flag=true 
		flag=MatchingRow(OriginalData,UniqueData,i,j,p,flag)
		if flag
			RowNos[i]=j
		end
	end	
end
RowNos=Int.(RowNos)
end