function EdgeConstraints(P1,P2,P3,MidPoint)
# EdgeCons: 
#               
# usage:
# [TriNo]=EdgeCons(Pa,Pb);
#
# Arguments: (input)
# Pa,Pb          - The corner point of each triangle in 'Triangles'.
#                    Arranged so the row index's correspond exactly to
#                    'Triangles' and 'MidPoint'. 
#
#
# Arguments: (output)
# TriNo           - The triangle index that is the connected triangle to
#                   that edge. Location of index is the triangle in question.
#
# P1P2FreeFlg     
# P2P3FreeFlg
# P3P1FreeFlg     - Flags of free edges between Pa and Pb (location
#                   specific)
#                             
# Example usage:
#
# 
#  Author: Tim Davis
#  Copyright 2017, Tim Davis, Potsdam University

#Checking 6 points
(SixPntsP1P2)=CreateSortedEdgePoints(P1,P2);
(SixPntsP2P3)=CreateSortedEdgePoints(P2,P3);
(SixPntsP3P1)=CreateSortedEdgePoints(P3,P1);

#Flags for free edges:
P1P2FreeFlg=fill(true,length(P1[:,1])); #zeros(length(P1[:,1]));
P2P3FreeFlg=fill(true,length(P1[:,1]));
P3P1FreeFlg=fill(true,length(P1[:,1]));

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
for j=1:length(P1[:,1])
	for i=1:length(P1[:,1])
        
    #Check to see if its worth continuing, here if the distance between the
    #midpoints is further than the sum of the two triangle perimeters(/2) we
    #skip. 
    if MidDist[i,j]>(HPerimP[i]+HPerimP[j])
        #disp('skipping')
        continue #Skip to next iteration
    end
    
    #Diff edge cons:
    if MatchingRow(SixPntsP1P2,SixPntsP2P3,i,j) 
        #all(SixPntsP1P2[i,:]==SixPntsP2P3[j,:])
        P1P2FreeFlg[i]=false; 
        P2P3FreeFlg[j]=false; 
    end
    if MatchingRow(SixPntsP1P2,SixPntsP3P1,i,j) 
        #all(SixPntsP1P2[i,:]==SixPntsP3P1[j,:])
        P1P2FreeFlg[i]=false; 
        P3P1FreeFlg[j]=false; 
    end
    if MatchingRow(SixPntsP2P3,SixPntsP3P1,i,j) 
        #all(SixPntsP2P3[i,:]==SixPntsP3P1[j,:])
        P2P3FreeFlg[i]=false;  
        P3P1FreeFlg[j]=false; 
    end
        
    if i==j
        continue
    end
    
    #Self cons: (e.g. only P1P2 to P1P2 edges)
    if MatchingRow(SixPntsP1P2,SixPntsP1P2,i,j) 
        #all(SixPntsP1P2[i,:]==SixPntsP1P2[j,:])
        P1P2FreeFlg[i]=false; 
        P1P2FreeFlg[j]=false; 
    end
    
    if MatchingRow(SixPntsP2P3,SixPntsP2P3,i,j) 
        #all(SixPntsP2P3[i,:]==SixPntsP2P3[j,:])
        P2P3FreeFlg[i]=false; 
        P2P3FreeFlg[j]=false; 
    end
    
    if MatchingRow(SixPntsP3P1,SixPntsP3P1,i,j) 
        #all(SixPntsP3P1[i,:]==SixPntsP3P1[j,:])
        P3P1FreeFlg[i]=false; 
        P3P1FreeFlg[j]=false; 
    end

	end
end

return P1P2FreeFlg,P2P3FreeFlg,P3P1FreeFlg

end



