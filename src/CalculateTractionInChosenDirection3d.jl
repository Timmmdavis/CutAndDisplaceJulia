function CalculateTractionInChosenDirection3d( Tx,Ty,Tz,CosAx,CosAy,CosAz,ChosenDirectionCos )
# CalculateTractionInChosenDirection3d: Calculates traction on a plane in 
#                   any direction the user chooses. The user just needs to
#                   supply the direction cosines of the direction they want
#                   this in and the Cartesian traction components on this
#                   plane. If needed these can be calculated from a full
#                   stress tensor on the plane using function:
#                   "TractionVectorCartesianComponents3d.m".
#                   Based on Equation 6.52 in Pollard, D.D. and Fletcher,
#                   R.C., 2005. Fundamentals of structural geology.
#                   Cambridge University Press.
#               
# usage #1:
# [ TractionInDirection ] = CalculateTractionInChosenDirection3d( Tx,Ty,Tz,CosAx,CosAy,CosAz,ChosenDirectionCos )
#
# Arguments: (input)
# Tx,Ty,Tz          - The Cartesian traction components on the plane.
#
# CosAx,CosAy,CosAz - The seperated direction cosines of the surface 
#                    (column vectors)
#
# ChosenDirectionCos- 3*n array of the direction the user desired the
#                    traction in. This is direction cosines
#                    [CosAx,CosAy,CosAz] of a vector that faces in any
#                    direction.
#
# Arguments: (output)
# TractionInDirection- Column vector of traction magnitudes in the chosen
#                      direction.
#
# Example usage 1:
#
# # Calc dip slip traction for a plane dipping 45 degrees facing east:
# CosAx=cosd(45);
# CosAy=0;
# CosAz=cosd(45);
# Sxx=1; Syy=0; Szz=0; Sxy=0; Sxz=0; Syz=0;
# [ Tx,Ty,Tz ] = TractionVectorCartesianComponents3d(  Sxx,Syy,Szz,Sxy,Sxz,Syz,CosAx,CosAy,CosAz );
# [ StrikeSlipCosine,DipSlipCosine ] = CalculateDSandSSDirs( [CosAx,CosAy,CosAz] ) ;  
# [ DipSlipTraction ] = CalculateTractionInChosenDirection3d( Tx,Ty,Tz,CosAx,CosAy,CosAz,DipSlipCosine )
#
#  Author: Tim Davis
#  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


#Init vector to fill
TractionInDirection=zeros(size(Tx));
for i=1:size(Tx,1) #On each col

	CosAx2=1-(CosAx[i].^2);
	CosAy2=1-(CosAy[i].^2);
	CosAz2=1-(CosAz[i].^2);
	CosAxy=CosAx[i]*CosAy[i];
	CosAyz=CosAy[i]*CosAz[i];
	CosAzx=CosAz[i]*CosAx[i];

	for j=1:size(Tx,2) #Running through the rows..

		#Split up the equation into seperate parts (each line in the book):
		A=(Tx[i,j]* CosAx2 - Ty[i,j]*CosAxy -Tz[i,j]*CosAzx)*ChosenDirectionCos[i,1];
		B=(Tx[i,j]*-CosAxy + Ty[i,j]*CosAy2 -Tz[i,j]*CosAyz)*ChosenDirectionCos[i,2];
		C=(Tx[i,j]*-CosAzx - Ty[i,j]*CosAyz +Tz[i,j]*CosAz2)*ChosenDirectionCos[i,3];
		#Sum the parts of the equation
		TractionInDirection[i,j]=A+B+C; 
		
	end
end	

return TractionInDirection 

end

