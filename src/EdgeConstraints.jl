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

println("checking 6pnts")
(SixPntsP1P2)=CreateSortedEdgeVec(P1,P2);
(SixPntsP2P3)=CreateSortedEdgeVec(P2,P3);
(SixPntsP3P1)=CreateSortedEdgeVec(P3,P1);
println("finished check")

#Flags for free edges:
P1P2FreeFlg=zeros(length(P1[:,1]));
P2P3FreeFlg=zeros(length(P1[:,1]));
P3P1FreeFlg=zeros(length(P1[:,1]));

#Calculate half length of each triangles perimeter
( ~,HPerimP ) = AreaOfTriangle3D( P1[:,1],P1[:,2],P1[:,3],P2[:,1],P2[:,2],P2[:,3],P3[:,1],P3[:,2],P3[:,3] );
#Calculate the distance between all midpoints:
num=length(MidPoint[:,1]);
X=repeat(MidPoint[:,1],1,num)-repeat(MidPoint[:,1]',num,1);
Y=repeat(MidPoint[:,2],1,num)-repeat(MidPoint[:,2]',num,1);
Z=repeat(MidPoint[:,3],1,num)-repeat(MidPoint[:,3]',num,1);
MidDist=sqrt.((X.^2)+(Y.^2)+(Z.^2));


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
    if all(SixPntsP1P2[i,:]==SixPntsP2P3[j,:])
        P1P2FreeFlg[i,:]=P1P2FreeFlg[i,:].+1; 
        P2P3FreeFlg[j,:]=P2P3FreeFlg[j,:].+1; 
    end
    
    if all(SixPntsP1P2[i,:]==SixPntsP3P1[j,:])
        P1P2FreeFlg[i,:]=P1P2FreeFlg[i,:].+1; 
        P3P1FreeFlg[j,:]=P3P1FreeFlg[j,:].+1; 
    end
    
    if all(SixPntsP2P3[i,:]==SixPntsP3P1[j,:])
        P2P3FreeFlg[i,:]=P2P3FreeFlg[i,:].+1;  
        P3P1FreeFlg[j,:]=P3P1FreeFlg[j,:].+1; 
    end
        
    if i==j
        continue
    end
    
    #Self cons: (e.g. only P1P2 to P1P2 edges)
    #if isequal(SixPntsP1P2[i,:],SixPntsP1P2[j,:])
    if all(SixPntsP1P2[i,:]==SixPntsP1P2[j,:])
        P1P2FreeFlg[i,:]=P1P2FreeFlg[i,:].+1; 
        P1P2FreeFlg[j,:]=P1P2FreeFlg[j,:].+1; 
    end
    if all(SixPntsP2P3[i,:]==SixPntsP2P3[j,:])
        P2P3FreeFlg[i,:]=P2P3FreeFlg[i,:].+1; 
        P2P3FreeFlg[j,:]=P2P3FreeFlg[j,:].+1; 
    end
    if all(SixPntsP3P1[i,:]==SixPntsP3P1[j,:])
        P3P1FreeFlg[i,:]=P3P1FreeFlg[i,:].+1; 
        P3P1FreeFlg[j,:]=P3P1FreeFlg[j,:].+1; 
    end

	end
end

#Making flags
P1P2FreeFlg=P1P2FreeFlg.==0;
P2P3FreeFlg=P2P3FreeFlg.==0;
P3P1FreeFlg=P3P1FreeFlg.==0;

return P1P2FreeFlg,P2P3FreeFlg,P3P1FreeFlg

end

function CreateSortedEdgeVec(Pa,Pb)
SixPnts=zeros(length(Pa[:,1]),6);
for i=1:length(Pa[:,1])

	#If the two are not equal in X
	if Pa[i,1]!=Pb[i,1]
		#X for Pa is bigger
        if Pa[i,1]<Pb[i,1]
			SixPnts[i,:]=[Pa[i,:] Pb[i,:]];
        else
			SixPnts[i,:]=[Pb[i,:] Pa[i,:]];
        end
	#So we look at Y	
	elseif Pa[i,2]!=Pb[i,2]
		#Y for Pa is bigger
        if Pa[i,2]<Pb[i,2]
			SixPnts[i,:]=[Pa[i,:] Pb[i,:]];
        else
			SixPnts[i,:]=[Pb[i,:] Pa[i,:]];
        end	
	#So we look at Z			
    elseif Pa[i,3]!=Pb[i,3]
		#Z for Pa is bigger
        if Pa[i,3]<Pb[i,3]
			SixPnts[i,:]=[Pa[i,:] Pb[i,:]];
        else
			SixPnts[i,:]=[Pb[i,:] Pa[i,:]];
        end		
	else
		#Points are equal anyway
	SixPnts[i,:]=[Pa[i,:] Pb[i,:]];
	end
	
end

return SixPnts

end
# #rounding if needed
# roundV=10;
# SixPnts=round(SixPnts,roundV);

