#Test case checking that the elastic constants functions returns the right values

GIn=ShearModulus(50); 
νIn=PoissonsRatio(0.01);

@info GIn.G
@info νIn.ν

(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(GIn,νIn);
K=BulkModulus(K);E=YoungsModulus(E);

(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(K,E);
K=BulkModulus(K);λ=LamesConstant(λ)

(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(K,λ);
K=BulkModulus(K);ν=PoissonsRatio(ν)

(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(K,ν);
K=BulkModulus(K);G=ShearModulus(G)

(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(K,G);
E=YoungsModulus(E);λ=LamesConstant(λ)

(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(E,λ);
E=YoungsModulus(E);ν=PoissonsRatio(ν)

(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(E,ν);
E=YoungsModulus(E);G=ShearModulus(G)

(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(E,G);
λ=LamesConstant(λ);ν=PoissonsRatio(ν)

(K,E,λ,ν,G) = CutAndDisplaceJulia.ElasticConstantsCheck(λ,ν);
λ=LamesConstant(λ);G=ShearModulus(G)

(K,E,λ,ν2,G2) = CutAndDisplaceJulia.ElasticConstantsCheck(λ,G);

@info GIn.G 	G2
@info νIn.ν	ν2

if !isapprox(GIn.G,G2) || !isapprox(νIn.ν,ν2)
	error("Error! Values not equal after passing through the elastic constants check")
end

println("Test Passed")

