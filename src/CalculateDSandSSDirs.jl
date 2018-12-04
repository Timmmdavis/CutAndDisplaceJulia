function CalculateDSandSSDirs(CosAx,CosAy,CosAz)
# CalculateDSandSSDirs: Calculates the dipslip and strikeslip direction
#                    cosines from a column normal vector with 3 cols
#                    representing aX aY aZ.
#                    Figure drawing can be turned on at the end of the
#                    function to see these directions. 
#
#               
#
# usage #1:
# [ StrikeSlipCosine,DipSlipCosine ] = CalculateDSandSSDirs( [],CosAx,CosAy,CosAz )
#
# Arguments: (input)
#
# CosAx,CosAy,CosAz - The seperated direction cosines of the surface 
#                    (column vectors)
#						CosAx=FaceNormalVector(:,1); 
#						CosAy=FaceNormalVector(:,2);
#						CosAz=FaceNormalVector(:,3);
#
# Arguments: (output)
# StrikeSlipCosine  - (n*3) array of CosAx,CosAy,CosAz directions of the 
#                     strike direction of the input plane
#
# DipSlipCosine     - (n*3) array of CosAx,CosAy,CosAz directions of the 
#                     dip direction of the input plane
#
# Example usage 1:
#
# # Get cosines for a plane dipping 45 degrees facing east:
# ( StrikeSlipCosine,DipSlipCosine ) = CalculateDSandSSDirs( cosd(45),0,cosd(45) )
#
#  Author: Tim Davis
#  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen


#Strikeslip
#Creating direction cosines for a vector pointing along the strike slip
#direction of each triangle. This is rotated so it represents right hand rule i.e.. strike
#is 90 deg to the right of the dip azimuth. 
#Just the 2D components of the normal vector (looking from above)
Res2D=[CosAx CosAy zeros(size(CosAx))];
#Normalise these (each row)
Res2D=normr(Res2D);

#Find vector that is 90 to this in 2D (XY)
StrikeSlipCosine=[-Res2D[2,1] Res2D[1,:] zeros(size(CosAx))];

#Dipslip
#Creating direction cosines for a vector pointing down the DipSlip
#direction of each triangle. 
#Calculate the new CosAz vector
DSCosAz=cos(acos(CosAz)-(pi/2));    #(90 deg to original)
#Calculate the length of vector NxNy when looked above in XY
L=sin(acos(DSCosAz));     
downdip=CosAz<0; #CosAz points down so we flip the length sign
if any(downdip)
	L[downdip]=-L[downdip];
end
#Calculate the angle a that NxNy vectors face in, this doesnt change when
#vector is rotated around dip  
a=atan(CosAy,CosAx);
#Calculate new lengths
DSCosAx=cos(a).*L;
DSCosAy=sin(a).*L;
DipSlipCosine=[-DSCosAx -DSCosAy DSCosAz];



# Now we make sure flat triangles follow the same conv as in Mehdi Nikko's
# TDE Script that is: "# Calculate unit strike, dip and normal to TD
# vectors: For a horizontal TD as an exception, if the normal vector points
# upward, the strike and dip vectors point Northward and Westward, whereas
# if the normal vector points downward, the strike and dip vectors point
# Southward and Westward, "

#Mehdi Nikkhoo's check for flat tris
eZ = [0;0;1];
for i=1:length(CosAx)
Vstrike = cross(eZ,[CosAx[i];CosAy[i];CosAz[i]]);
    #Applying a convention if flat: 
    if norm(Vstrike)==0
        if CosAz[i,:]<1 #if the surface norm faces up
            DipSlipCosine[i,:]   = [-1 0 0]; #pointing west
            StrikeSlipCosine[i,:]= [0 -1 0]; #pointing south 
        else
            DipSlipCosine[i,:]   = [-1 0 0]; #pointing west
            StrikeSlipCosine[i,:]= [0 1 0];  #pointing north  
        end

    end
end


return(StrikeSlipCosine,DipSlipCosine)
end
