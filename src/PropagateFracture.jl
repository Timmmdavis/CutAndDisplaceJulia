function PropagateFracture( FeP1P2S,FeP1P3S,FeP2P3S,FaceNormalVector,G,ν,KCrit,AvgTriangleEdgeLength )
#PropagateFracture 

#Find index locations:
IndxP1P2=findall(FeP1P2S.FreeFlg);
IndxP1P3=findall(FeP1P3S.FreeFlg);
IndxP2P3=findall(FeP2P3S.FreeFlg);

#Find MaxStrainEnergy value for all edges
MaxStrainEnergy=maximum(filter(!isnan,[FeP1P2S.StrainEnergy; FeP1P3S.StrainEnergy; FeP2P3S.StrainEnergy]))

(PV1,Ang1)=CreateNewEdgePoint(FeP1P2S,IndxP1P2,FaceNormalVector,G,ν,KCrit,MaxStrainEnergy,AvgTriangleEdgeLength)
(PV2,Ang2)=CreateNewEdgePoint(FeP1P3S,IndxP1P3,FaceNormalVector,G,ν,KCrit,MaxStrainEnergy,AvgTriangleEdgeLength)
(PV3,Ang3)=CreateNewEdgePoint(FeP2P3S,IndxP2P3,FaceNormalVector,G,ν,KCrit,MaxStrainEnergy,AvgTriangleEdgeLength)



return PV1,PV2,PV3,Ang1,Ang2,Ang3
end

function CreateNewEdgePoint(Fe,Indx,FaceNormalVector,G,ν,KCrit,MaxStrainEnergy,AvgTriangleEdgeLength)
	PointVectorX=[];
    PointVectorY=[];
    PointVectorZ=[];
    AngVector=[];
	for i=1:length(Indx)

        
	    #Extract some values
        I=Indx[i]; 
        #Check if the crack tip will extend
        if Fe.StrainEnergy[I]>KCrit #Plane strain criteria
            ( NwPntCX,NwPntCY,NwPntCZ,Ang ) = FindPropAngleAndPoint( Fe,FaceNormalVector,I,MaxStrainEnergy,AvgTriangleEdgeLength );
        else 
            continue
        end
        PointVectorX=push!(PointVectorX,NwPntCX)
        PointVectorY=push!(PointVectorY,NwPntCY)
        PointVectorZ=push!(PointVectorZ,NwPntCZ)
        AngVector=push!(AngVector,Ang)
    end
    PointVector=[PointVectorX PointVectorY PointVectorZ]
    return PointVector,AngVector
end

