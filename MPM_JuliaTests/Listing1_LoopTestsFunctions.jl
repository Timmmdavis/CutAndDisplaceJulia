m=1
n=10000 #000

#create and assign arrays
#values of m and n defined in Table 1

A01 = zeros(Float64, m, n);
A02 = zeros(Float64, m, n);
A03 = zeros(Float64, m, n);

#A01 = Array{Float64}(m, n) # m by n 2D array of Float64 numbers
#A02 = Array{Float64}(undef, m, n)
#A03 = Array{Float64}(undef, m, n)
for i = 1:1:m
	for j = 1:1:n
		A01[i,j] = 1.0*i
		A02[i,j] = -1.0*i
		A03[i,j] = 0
	end
end

Iterations = 1000


# using for i-j loop
tic=time()
for count in 1:1:Iterations
	Loop_ij(A01,A02,A03,m,n)
end
toc=time()
println("Elapsed time")
println(Iterations/(toc-tic))

 # using for j-i loop , should be used
tic=time()
for count in 1:1:Iterations
	for j=1:1:n
		for i=1:1:m
			A03[i,j] = A01[i,j]+A02[i,j]
		end
	end
end
toc=time()
println("Elapsed time")
println(Iterations/(toc-tic))


 # using vectorization
tic=time()
for count in 1:1:Iterations
	A03 = A01+A02
end
toc=time()
println("Elapsed time")
println(Iterations/(toc-tic))


function Loop_ij(A01,A02,A03,m,n)
	for i=1:1:m
		for j=1:1:n
			A03[i,j] = A01[i,j]+A02[i,j]
		end
	end
end