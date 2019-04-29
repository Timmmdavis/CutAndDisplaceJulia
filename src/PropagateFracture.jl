function PropagateFracture( FeP1P2S,FeP1P3S,FeP2P3S,FaceNormalVector,G,ν,KCrit )
#PropagateFracture 

#Find index locations:
IndxP1P2=findall(FeP1P2S.FreeFlg);
IndxP1P3=findall(FeP1P3S.FreeFlg);
IndxP2P3=findall(FeP2P3S.FreeFlg);

(PV1,Ang1)=CreateNewEdgePoint(FeP1P2S,IndxP1P2,FaceNormalVector,G,ν,KCrit)
(PV2,Ang2)=CreateNewEdgePoint(FeP1P3S,IndxP1P3,FaceNormalVector,G,ν,KCrit)
(PV3,Ang3)=CreateNewEdgePoint(FeP2P3S,IndxP2P3,FaceNormalVector,G,ν,KCrit)
return PV1,PV2,PV3,Ang1,Ang2,Ang3
end

function CreateNewEdgePoint(Fe,Indx,FaceNormalVector,G,ν,KCrit)
	PointVectorX=[];
    PointVectorY=[];
    PointVectorZ=[];
    AngVector=[];
	for i=1:length(Indx)
        #=
        if i==4
            @bp
                #Extract some values
            I=Indx[i]; 
            NrmVec=FaceNormalVector[I,:];
            #Check if the crack tip will extend
            ( StrainEnergy ) = StrainEnergyRelease(Fe.K1[I],Fe.K2[I],Fe.K3[I],G,ν);
            #@info StrainEnergy
            @enter FindPropAngleAndPoint( Fe.FeMd[I,:],Fe.FeM2Ev[I,:],Fe.FeLe[I],Fe.FeEv[I,:],NrmVec,Fe.K2[I],Fe.K1[I] );
        end
        =#
        
	    #Extract some values
        I=Indx[i]; 
        NrmVec=FaceNormalVector[I,:];
        #Check if the crack tip will extend
        ( StrainEnergy ) = StrainEnergyRelease(Fe.K1[I],Fe.K2[I],Fe.K3[I],G,ν);
        #@info StrainEnergy

        if StrainEnergy>KCrit #Plane strain criteria
            ( NwPntCX,NwPntCY,NwPntCZ,Ang ) = FindPropAngleAndPoint( Fe.FeMd[I,:],Fe.FeM2Ev[I,:],
                Fe.FeLe[I],Fe.FeEv[I,:],NrmVec,Fe.K2[I],Fe.K1[I] );
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

