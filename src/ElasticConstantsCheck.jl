#function ElasticConstantsCheck()
#ElasticConstantsCheck: Function to check the elastic constants input into
#               the function are consistent. Then allows for export of
#               constants that may be required. 
#               Input arguments must be predefined with the correct name
#               (see inputs description). This is case insensitive, the
#               inputs must also lie within reasonable elastic bounds.
# 
#
# usage #1.0: Get all elastic constants from two input
#Example
#G=ShearModulus(2)
#ν=PoissonsRatio(0.25)
#In any arg order
#(K,E,λ,ν,G)=CutAndDisplaceJulia.ElasticConstantsCheck(G=2),ν=ν)
#
#               
#
# Arguments: (input)
#   K           - Bulk modulus
#
#   E           - Young's modulus
#
# 	λ         - Lame's constant
#
#   ν         - Poisson's ratio 
#
#   G         - Shear modulus
#
#
#
#  Author: Tim Davis
#  Copyright 2019 Tim Davis, Potsdam University



function ElasticConstantsCheck(G::ShearModulus,ν::PoissonsRatio)
    
	#Extract out of structs (multiple dispatch)
	G=G.G;
	ν=ν.ν;
	#Equations
    K=((2*G)*(1+ν))/(3*(1-(2*ν)));
    E=(2*G)*(1+ν);
    λ=(2*G*ν)/(1-(2*ν));
	#Return results
	return (K,E,λ,ν,G)
	
end

function ElasticConstantsCheck(K::BulkModulus,E::YoungsModulus)
    
	#Extract out of structs (multiple dispatch)
	K=K.K;
	E=E.E;
	#Equations
    λ=((3*K)*((3*K)-E))/((9*K)-E);
    ν=((3*K)-E)/(6*K);
    G=(3*K*E)/((9*K)-E);
	#Return results
	return (K,E,λ,ν,G)
	
end

function ElasticConstantsCheck(K::BulkModulus,λ::LamesConstant)

	#Extract out of structs (multiple dispatch)
    K=K.K;
	λ=λ.λ;
	#Equations
    E=((9*K)*(K-λ))/((3*K)-λ);
    ν=λ/((3*K)-λ);
    G=(3/2)*(K-λ);
	#Return results	
	return (K,E,λ,ν,G)	
	
end

function ElasticConstantsCheck(K::BulkModulus,ν::PoissonsRatio)
	
	#Extract out of structs (multiple dispatch)
    K=K.K;
	ν=ν.ν;
	#Equations
    E=(3*K)*(1-(2*ν));
    λ=(3*K*ν)/(1+ν);
    G=((3*K)*(1-(2*ν)))/(2*(1+ν));
	#Return results
	return (K,E,λ,ν,G)		
	
end

function ElasticConstantsCheck(K::BulkModulus,G::ShearModulus)

	#Extract out of structs (multiple dispatch)
    K=K.K;
	G=G.G;
	#Equations
    E=(9*K*G)/((3*K)+G);
    λ=((3*K)-(2*G))/3; 
    ν=((3*K)-(2*G))/((6*K)+(2*G)); 
	#Return results
	return (K,E,λ,ν,G)		
	
end

function ElasticConstantsCheck(E::YoungsModulus,λ::LamesConstant)

	#Extract out of structs (multiple dispatch)
    E=E.E;
	λ=λ.λ;
	#Equations
    R=sqrt((E^2)+(9*(λ^2))+(2*E*λ));
    K=(E+(3*λ)+R)/6;
    ν=(2*λ)/(E+λ+R);
    G=(E-(3*λ)+R)/4;
	#Return results
	return (K,E,λ,ν,G)		
	
end

function ElasticConstantsCheck(E::YoungsModulus,ν::PoissonsRatio)

	#Extract out of structs (multiple dispatch)
    E=E.E;
	ν=ν.ν;
	#Equations
    K=E/(3*(1-(2*ν)));
    λ=(E*ν)/((1+ν)*(1-(2*ν)));
    G=E/(2*(1+ν));
	#return results
	return (K,E,λ,ν,G)		
	
end

function ElasticConstantsCheck(E::YoungsModulus,G::ShearModulus)

	#Extract out of structs (multiple dispatch)
    E=E.E;
	G=G.G;
	#Equations
    K=(E*G)/(3*((3*G)-E));
    λ=(G*(E-(2*G)))/((3*G)-E);
    ν=(E-(2*G))/(2*G);
	#return results
	return (K,E,λ,ν,G)		
	
end

function ElasticConstantsCheck(λ::LamesConstant,ν::PoissonsRatio)

	#Extract out of structs (multiple dispatch)
    λ=λ.λ
	ν=ν.ν;
	#Equations
    K=(λ*(1+ν))/(3*ν);
    E=(λ*(1+ν)*(1-(2*ν)))/ν;
    G=(λ*(1-(2*ν)))/(2*ν);
	#return results
	return (K,E,λ,ν,G)	 
	
end

function ElasticConstantsCheck(λ::LamesConstant,G::ShearModulus)

	#Extract out of structs (multiple dispatch)
    λ=λ.λ
	G=G.G;
	#Equations
    K=((3*λ)+(2*G))/3;
    E=(G*((3*λ)+(2*G)))/(λ+G);
    ν=λ/(2*(λ+G));
	#return results
	return (K,E,λ,ν,G)		
	
end

######################################################

#Defining functions so input order is independant
function ElasticConstantsCheck(ν::PoissonsRatio,G::ShearModulus)
	(K,E,λ,ν,G)=ElasticConstantsCheck(G::ShearModulus,ν::PoissonsRatio)
	return (K,E,λ,ν,G)
end

function ElasticConstantsCheck(E::YoungsModulus,K::BulkModulus,)
	(K,E,λ,ν,G)=ElasticConstantsCheck(K::BulkModulus,E::YoungsModulus)
	return (K,E,λ,ν,G)
end

function ElasticConstantsCheck(λ::LamesConstant,K::BulkModulus)
	(K,E,λ,ν,G)=ElasticConstantsCheck(K::BulkModulus,λ::LamesConstant)
	return (K,E,λ,ν,G)	
end

function ElasticConstantsCheck(ν::PoissonsRatio,K::BulkModulus,)
	(K,E,λ,ν,G)=ElasticConstantsCheck(K::BulkModulus,ν::PoissonsRatio)
	return (K,E,λ,ν,G)		
end

function ElasticConstantsCheck(G::ShearModulus,K::BulkModulus,)
	(K,E,λ,ν,G)=ElasticConstantsCheck(K::BulkModulus,G::ShearModulus)
	return (K,E,λ,ν,G)		
end

function ElasticConstantsCheck(λ::LamesConstant,E::YoungsModulus)
	(K,E,λ,ν,G)=ElasticConstantsCheck(E::YoungsModulus,λ::LamesConstant)
	return (K,E,λ,ν,G)		
end

function ElasticConstantsCheck(ν::PoissonsRatio,E::YoungsModulus)
	(K,E,λ,ν,G)=ElasticConstantsCheck(E::YoungsModulus,ν::PoissonsRatio)
	return (K,E,λ,ν,G)		
end

function ElasticConstantsCheck(G::ShearModulus,E::YoungsModulus)
	(K,E,λ,ν,G)=ElasticConstantsCheck(E::YoungsModulus,G::ShearModulus)
	return (K,E,λ,ν,G)		
end

function ElasticConstantsCheck(ν::PoissonsRatio,λ::LamesConstant)
	(K,E,λ,ν,G)=ElasticConstantsCheck(λ::LamesConstant,ν::PoissonsRatio)
	return (K,E,λ,ν,G)
end

function ElasticConstantsCheck(G::ShearModulus,λ::LamesConstant)
	(K,E,λ,ν,G)=ElasticConstantsCheck(λ::LamesConstant,G::ShearModulus)
	return (K,E,λ,ν,G)		
end