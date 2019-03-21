function normr(A) 
#using LinearAlgebra
#inputs like:    Res2D=[CosAx CosAy zeros(size(CosAx))];

B= Array{Float64}(undef, size(A)); 
if size(A,1)==1
	B=normalize(vec(A))
else
	
	for i=1:size(A,1)
		B[i,:]=normalize(A[i,:],2)
	end
end

return(B)
end