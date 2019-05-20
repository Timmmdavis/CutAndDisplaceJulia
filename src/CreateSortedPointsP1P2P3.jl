function CreateSortedPointsP1P2P3(Pa,Pb,Pc)
#Used like:
##Checking 6 points
#(SixPntsP1P2)=CreateSortedEdgePoints(P1,P2);
#(SixPntsP2P3)=CreateSortedEdgePoints(P2,P3);
#(SixPntsP3P1)=CreateSortedEdgePoints(P3,P1);


NinePnts=zeros(length(Pa[:,1]),9);
for i=1:length(Pa[:,1])

	#If the two are not equal in X
	if Pa[i,1]!=Pb[i,1]!=Pc[i,1]
		#X for Pa is bigger
        if Pa[i,1]<Pb[i,1]
            if Pa[i,1]<Pc[i,1]
                if Pb[i,1]<Pc[i,1]
                    #Pa < Pb < Pc
                    NinePnts=AddToVect2(NinePnts,Pa,Pb,Pc,i); 
                else 
                    #Pa < Pc < Pb
                    NinePnts=AddToVect2(NinePnts,Pa,Pc,Pb,i); 
                end
            else #Pc<Pa
                #Pc < Pa < Pb
                NinePnts=AddToVect2(NinePnts,Pc,Pa,Pb,i); 
            end
        else #Pb<Pa
            if Pb[i,1]<Pc[i,1]
                if Pa[i,1]<Pc[i,1]
                    NinePnts=AddToVect2(NinePnts,Pb,Pa,Pc,i); 
                else 
                    NinePnts=AddToVect2(NinePnts,Pb,Pc,Pa,i); 
                end
            else #Pc<Pb & Pb<Pa
                 NinePnts=AddToVect2(NinePnts,Pc,Pb,Pa,i); 
            end
        end
	#So we look at Y	
    elseif Pa[i,2]!=Pb[i,2]!=Pc[i,2]
        #Y for Pa is bigger
        if Pa[i,2]<Pb[i,2]
            if Pa[i,2]<Pc[i,2]
                if Pb[i,2]<Pc[i,2]
                    #Pa < Pb < Pc
                    NinePnts=AddToVect2(NinePnts,Pa,Pb,Pc,i); 
                else 
                    #Pa < Pc < Pb
                    NinePnts=AddToVect2(NinePnts,Pa,Pc,Pb,i); 
                end
            else #Pc<Pa
                #Pc < Pa < Pb
                NinePnts=AddToVect2(NinePnts,Pc,Pa,Pb,i); 
            end
        else #Pb<Pa
            if Pb[i,2]<Pc[i,2]
                if Pa[i,2]<Pc[i,2]
                    NinePnts=AddToVect2(NinePnts,Pb,Pa,Pc,i); 
                else 
                    NinePnts=AddToVect2(NinePnts,Pb,Pc,Pa,i); 
                end
            else #Pc<Pb & Pb<Pa
                 NinePnts=AddToVect2(NinePnts,Pc,Pb,Pa,i); 
            end
        end
	#So we look at Z			
    elseif Pa[i,3]!=Pb[i,3]!=Pc[i,3]
        #Z for Pa is bigger
        if Pa[i,3]<Pb[i,3]
            if Pa[i,3]<Pc[i,3]
                if Pb[i,3]<Pc[i,3]
                    #Pa < Pb < Pc
                    NinePnts=AddToVect2(NinePnts,Pa,Pb,Pc,i); 
                else 
                    #Pa < Pc < Pb
                    NinePnts=AddToVect2(NinePnts,Pa,Pc,Pb,i); 
                end
            else #Pc<Pa
                #Pc < Pa < Pb
                NinePnts=AddToVect2(NinePnts,Pc,Pa,Pb,i); 
            end
        else #Pb<Pa
            if Pb[i,3]<Pc[i,3]
                if Pa[i,3]<Pc[i,3]
                    NinePnts=AddToVect2(NinePnts,Pb,Pa,Pc,i); 
                else 
                    NinePnts=AddToVect2(NinePnts,Pb,Pc,Pa,i); 
                end
            else #Pc<Pb & Pb<Pa
                 NinePnts=AddToVect2(NinePnts,Pc,Pb,Pa,i); 
            end
        end

    else
    println("Inside")
	#Points are equal anyway
    NinePnts=AddToVect2(NinePnts,Pa,Pb,Pc,i);

	end
	
end

return NinePnts

end

function AddToVect2(NinePnts,P1,P2,P3,i)
    for j=1:3
        NinePnts[i,j]=  P1[i,j];
        NinePnts[i,j+3]=P2[i,j];
        NinePnts[i,j+6]=P3[i,j];
    end
    #println(NinePnts[i,:])
    return NinePnts
end
# #rounding if needed
# roundV=10;
# SixPnts=round(SixPnts,roundV);