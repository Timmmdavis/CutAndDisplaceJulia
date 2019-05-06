function CreateSortedEdgePoints(Pa,Pb)
#Used like:
##Checking 6 points
#(SixPntsP1P2)=CreateSortedEdgePoints(P1,P2);
#(SixPntsP2P3)=CreateSortedEdgePoints(P2,P3);
#(SixPntsP3P1)=CreateSortedEdgePoints(P3,P1);


SixPnts=zeros(length(Pa[:,1]),6);
for i=1:length(Pa[:,1])

	#If the two are not equal in X
	if Pa[i,1]!=Pb[i,1]
		#X for Pa is bigger
        if Pa[i,1]<Pb[i,1]
			#SixPnts[i,:]=[Pa[i,:] Pb[i,:]];
            SixPnts=AddToVect(SixPnts,Pa,Pb,i);            
        else
			#SixPnts[i,:]=[Pb[i,:] Pa[i,:]];
            SixPnts=AddToVect(SixPnts,Pb,Pa,i);         
        end
	#So we look at Y	
	elseif Pa[i,2]!=Pb[i,2]
		#Y for Pa is bigger
        if Pa[i,2]<Pb[i,2]
			#SixPnts[i,:]=[Pa[i,:] Pb[i,:]];
            SixPnts=AddToVect(SixPnts,Pa,Pb,i);
        else
            #SixPnts[i,:]=[Pb[i,:] Pa[i,:]];
            SixPnts=AddToVect(SixPnts,Pb,Pa,i);
        end	
	#So we look at Z			
    elseif Pa[i,3]!=Pb[i,3]
		#Z for Pa is bigger
        if Pa[i,3]<Pb[i,3]
			#SixPnts[i,:]=[Pa[i,:] Pb[i,:]];
            SixPnts=AddToVect(SixPnts,Pa,Pb,i);
        else
            #SixPnts[i,:]=[Pb[i,:] Pa[i,:]];
            SixPnts=AddToVect(SixPnts,Pb,Pa,i);  
        end		
	else
	#Points are equal anyway
    #SixPnts[i,:]=[Pa[i,:] Pb[i,:]];
    SixPnts=AddToVect(SixPnts,Pa,Pb,i);

	end
	
end

return SixPnts

end

function AddToVect(SixPnts,P1,P2,i)
    for j=1:3
        SixPnts[i,j]=P1[i,j];
        SixPnts[i,j+3]=P2[i,j];
    end
    return SixPnts
end
# #rounding if needed
# roundV=10;
# SixPnts=round(SixPnts,roundV);