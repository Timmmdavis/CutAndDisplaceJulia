function ConnectedConstraints(P1,P2,P3,MidPoint)


#Checking 6 points
(SixPntsP1P2)=CreateSortedEdgePoints(P1,P2);
(SixPntsP2P3)=CreateSortedEdgePoints(P2,P3);
(SixPntsP3P1)=CreateSortedEdgePoints(P3,P1);

#rows nos are index of tri, the three rows are filled with index's of
#connected tris
#SortedTris=Int.(zeros(length(P1[:,1]),3)); #zeros(length(P1[:,1]));
SortedTris=fill(0,length(P1[:,1]),3)
#Col1=P1P2 Col2=P2P3 Col3=P1P3 - numbers are the connected tris edge row number
#i.e. 23 is edge P2P3 of the connected tri 
ConnectedEdge=fill(0,length(P1[:,1]),3)

#Calculate half length of each triangles perimeter
( ~,HPerimP ) = AreaOfTriangle3D( P1[:,1],P1[:,2],P1[:,3],P2[:,1],P2[:,2],P2[:,3],P3[:,1],P3[:,2],P3[:,3] );
#Calculate the distance between all midpoints:
num=length(MidPoint[:,1]);
X=repeat(MidPoint[:,1],1,num)-repeat(MidPoint[:,1]',num,1);
Y=repeat(MidPoint[:,2],1,num)-repeat(MidPoint[:,2]',num,1);
Z=repeat(MidPoint[:,3],1,num)-repeat(MidPoint[:,3]',num,1);
#Matrix of distances
MidDist=zeros(size(X))
for i=1:size(MidDist,1)
    for j=1:size(MidDist,2)
        MidDist[i,j]=sqrt((X[i,j]^2)+(Y[i,j]^2)+(Z[i,j]^2));
    end
end

#Do Pa Pb connections
for i=1:length(P1[:,1])
	for j=1:length(P1[:,1])
        
    #Check to see if its worth continuing, here if the distance between the
    #midpoints is further than the sum of the two triangle perimeters(/2) we
    #skip. 
    if MidDist[i,j]>(HPerimP[i]+HPerimP[j])
        #disp('skipping')
        continue #Skip to next iteration
    end
    
    #Diff edge cons:
    if MatchingRow(SixPntsP1P2,SixPntsP2P3,i,j) 
        #FillSortedTris(SortedTris,i,j)
        SortedTris[i,1]=j
        ConnectedEdge[i,1]=23;
    end
    if MatchingRow(SixPntsP1P2,SixPntsP3P1,i,j) 
        #FillSortedTris(SortedTris,i,j)
        SortedTris[i,1]=j
        ConnectedEdge[i,1]=13;
    end
    if MatchingRow(SixPntsP2P3,SixPntsP1P2,i,j) 
        #FillSortedTris(SortedTris,i,j)
        SortedTris[i,2]=j
        ConnectedEdge[i,2]=12
    end
    if MatchingRow(SixPntsP2P3,SixPntsP3P1,i,j) 
        #FillSortedTris(SortedTris,i,j)
        SortedTris[i,2]=j
        ConnectedEdge[i,2]=13
    end
    if MatchingRow(SixPntsP3P1,SixPntsP1P2,i,j) 
        #FillSortedTris(SortedTris,i,j)
        SortedTris[i,3]=j
        ConnectedEdge[i,3]=12
    end
    if MatchingRow(SixPntsP3P1,SixPntsP2P3,i,j) 
        #FillSortedTris(SortedTris,i,j)
        SortedTris[i,3]=j
        ConnectedEdge[i,3]=23
    end
        
    if i==j
        continue
    end
    
    #Self cons: (e.g. only P1P2 to P1P2 edges)
    if MatchingRow(SixPntsP1P2,SixPntsP1P2,i,j) 
        #FillSortedTris(SortedTris,i,j)
        SortedTris[i,1]=j
        ConnectedEdge[i,1]=12
    end
    
    if MatchingRow(SixPntsP2P3,SixPntsP2P3,i,j) 
        #FillSortedTris(SortedTris,i,j)
        SortedTris[i,2]=j
        ConnectedEdge[i,2]=23
    end
    
    if MatchingRow(SixPntsP3P1,SixPntsP3P1,i,j) 
        #FillSortedTris(SortedTris,i,j)
        SortedTris[i,3]=j
        ConnectedEdge[i,3]=13
    end

	end
end

return SortedTris,ConnectedEdge

end

function FillSortedTris(SortedTris,RowIndx,TriIndx)
    #RowIndx - the current tri 
    #TriIndx - known connection
    if SortedTris[RowIndx,1]==0
        SortedTris[RowIndx,1]=TriIndx;
    elseif SortedTris[RowIndx,2]==0
        SortedTris[RowIndx,2]=TriIndx;
    elseif SortedTris[RowIndx,3]==0
        SortedTris[RowIndx,3]=TriIndx; 
    end 
end      

function MatchingRow(A1,A2,i,j)
#Compare Els of two vectors
#memoryless version of: all(A1[i,:]==A2[j,:])
for p=1:6
    if A1[i,p]!=A2[j,p]
        return false
    end
end
return true
end



