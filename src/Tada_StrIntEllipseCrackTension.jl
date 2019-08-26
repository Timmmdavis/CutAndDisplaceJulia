function Tada_StrIntEllipseCrackTension(Szz,Radius_a,Radius_b,xLocCrackWall,yLocCrackWall)
# Tada_StrIntEllipseCrackTension: Returns the stress intensity 
#               factors at the tips of a elliptical shaped crack under tension.
#				Equation from: Tada Stress analysis handbook P.26.2:
#
# usage #1: 
# K1 = Tada_StrIntEllipseCrackTension(Szz,theta,Radius_a,Radius_b)
#
# Arguments: (input)
# Szz             - The remote uniaxial tension along the Z-axis.
#

#
# Radius_a        - The semiaxis radius a and b of the crack.
#
# Arguments: (output)
# K1              - K1 at the crack tips.
#
#
#
#  Author: Tim Davis
#  Copyright 2019, Tim Davis, Potsdam University


#using Elliptic
#Elliptic table (Tada: 654.pdf) : alphadeg=90=1 45=0.5 0=0. Get the same results for K and E
#E.g. Elliptic.E(0.5)

#=
Szz=1;
a=1; #A must be bigger than b
b=1; #Keep at 1
theta=1:1:360; 
theta=deg2rad.(theta)
=#




#Ellipse shaped crack: Tada etc P385(405.pdf) 26.2
a=Radius_a;
b=Radius_b;

# Theta           - The parametric angle on the crack wall, measured counter-clockwise 
#                   away from semiaxis a (in crack plane), radians. See Fig one page P.26.1
H=sqrt.((a^2).-(xLocCrackWall.^2) )
θ = atan.(H,xLocCrackWall)
for i=1:length(θ) 
	if yLocCrackWall[i]<0
		θ[i]=-θ[i]
	end
end

l=b.*sqrt.(sin.(θ).^2 .+((b^2)/(a^2)).*cos.(θ).^2); #Length to tip from centre (P384)
#l2=b.*((sin.(theta).^2 .+((b^2)/(a^2)).*cos.(theta).^2).^(1/4));

kSqrd=1-(b^2)/(a^2);

#Kk=Elliptic.K(k);
Ek=Elliptic.E(kSqrd);

K1=((Szz.*sqrt.(pi.*l))./Ek);


#=
#Should be the same when a and b are the same
K13Dpenny=(2/pi).*(Szz.*sqrt.(pi.*b))
#Should be the same when a=inf and theta=90/270
K12DLineCrack=(Szz.*sqrt.(pi.*b))

@info K1 K1simple K12DLineCrack
=#


return K1,θ

end