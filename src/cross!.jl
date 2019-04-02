#Cross from Julia Base.LinearAlgebra but with no allocations. 
function cross!(a::AbstractVector, b::AbstractVector,c::AbstractVector)
     
	#If you are using julias cross function put a and b to vecs first. i.e. a=vec(a)

     if !(length(a) == length(b) == length(c) == 3)
         throw(DimensionMismatch("cross product is only defined for vectors of length 3"))
     end
    c[1]=a[2]*b[3]-a[3]*b[2]
	c[2]=a[3]*b[1]-a[1]*b[3]
	c[3]=a[1]*b[2]-a[2]*b[1]


end