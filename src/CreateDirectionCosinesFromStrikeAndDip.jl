function CreateDirectionCosinesFromStrikeAndDip(Strike,Dip)
#Given strike and dip in degrees create direction cosines of the fault (normal,dip and strike). 

( StrikeSlipCosine,DipSlipCosine ) = CalculateDSandSSDirs( cosd(Dip),0,cosd(Dip) ) #Plane dipping at dip (facing east!)
RAngle=deg2rad(-Strike)
StrikeSlipCosine=RotateCosine3D(StrikeSlipCosine,RAngle,"z")
DipSlipCosine=RotateCosine3D(DipSlipCosine,RAngle,"z")
FaceNormalVector=cross(vec(StrikeSlipCosine),vec(DipSlipCosine))
FaceNormalVector=FaceNormalVector'
CheckDirectionCosinePerp(FaceNormalVector,DipSlipCosine,StrikeSlipCosine)

return(StrikeSlipCosine,DipSlipCosine,FaceNormalVector)

end