function normr(A) 
#using LinearAlgebra
#inputs like:    Res2D=[CosAx CosAy zeros(size(CosAx))];

#@info A

#B=normalize(vec(A))

B= Array{Float64}(undef, size(A)); 
if size(A,1)==1 
	B=LinearAlgebra.normalize(vec(A))
else
	
	for i=1:size(A,1)
		B[i,:]=LinearAlgebra.normalize(A[i,:],2)
	end
end

return(B)
end