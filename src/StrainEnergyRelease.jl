#StrainEnergyRelease: Calculates strain energy released 'Gc'. For use to
#                     see if this is over the critical fracture toughness.
#                     See Martel 2016: 'Progress in understanding sheeting
#                     joints over the past two centuries'. Eq.6
#               
# usage:
# [ StrainEnergy ] = StrainEnergyRelease(K1,K2,K3,G,ν)
#
# Arguments: (input)
# K1,K2,K3          - The stress intensity factors at the crack tip/edge. 
#
# G                - The shear modulus of the material
#
# ν                - The Young's modulus of the material
#
# Arguments: (output)
# StrainEnergy      - A value of the strain energy release rate. 
# fracturetoughness - A value that can be compared to Kcrit
#                             
# Example usage:
#
# N/a
# 
#  Author: Tim Davis
#  Copyright 2018, Tim Davis, Potsdam University

function StrainEnergyRelease(K1,K2,K3,G,ν)


#Eq 48 - Isogeometric boundary element methods for three dimensional static
#fracture and fatigue crack growth X. Peng, E. Atroshchenko,∗, P. Kerfriden, S.P.A. Bordas 2017
#Wiki G-Criterion - stress intensity factor page
E=(2.0*G)*(1.0+ν);
Eprime=E/(1-ν^2)
fracturetoughness=sqrt.((abs.(K1).^2)+(abs.(K2).^2)+((Eprime/(2*G))*(abs.(K3).^2)));


#=
#Tada G criterion p1.7 - Crack extension force
E=(2.0*G)*(1.0+ν);
Eprime=E/(1-ν^2)
G1=(abs.(K1).^2)./Eprime
G2=(abs.(K2).^2)./Eprime
G3=(abs.(K3).^2)./2*G
StrainEnergy=G1+G2+G3
=#

##Previous (same as wiki but simple )
#StrainEnergy=((abs.(K1).^2+abs.(K2).^2).*(1-ν)+abs.(K3).^2)./(2*G);

##From wikipedia strain erergy (same result)
#E=(2.0*G)*(1.0+ν);
#c1=(1.0-(ν^2.0))/E
#c2=1.0/(2.0*G)
#StrainEnergy=(abs.(K1).^2).*c1.+(abs.(K2).^2).*c1.+(abs.(K3).^2).*c2

##from wiki - rearranged to get fracturetoughness from G
#fracturetoughness=sqrt.(StrainEnergy.*Eprime);

return fracturetoughness
end

