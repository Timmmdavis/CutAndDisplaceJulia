#declarationofauserdefinedtype,named'DataStructure'
#withonemembernamedMwhichis1-dimensionalarray
#inwhatfollows,valuesofnandmvariedfollowingTable2
type DataStructure
	M::Arrayf{Float64} #thesymbol::usedtodefineaTYPE
	function DataStructure()
		new(zeros(n))
	end
end

#createandassigntypemembers
S01=Array{DataStructure}(m)
S02=Array{DataStructure}(m)
S03=Array{DataStructure}(m)
#thefollowingcanbesimplied:fori=1:m
for i=1:1:m
	S01[i]=DataStructure()
	S02[i]=DataStructure()
	S03[i]=DataStructure()
	for j=1:1:n
		S01[i].M[j]=1.0*i
		S02[i].M[j]=-1.0*i
		S03[i].M[j]=0
	end
end

#method1,usingcompositetypeobjectdirectly
tic=time()
for count in 1:1:iIterations
	fori=1:1:m
		forj=1:1:n
			S03[i].M[j]=S01[i].M[j]+S02[i].M[j]
		end
	end
end
toc=time()

#method2,usingreferencetocompositetypeobjects
tic=time()
for count in 1:1:iIterations
	for i=1:1:m
		this01=S01[i]
		this02=S02[i]
		this03=S03[i]
		for j=1:1:n
			this03.M[j]=this01.M[j]+this02.M[j]
		end
	end
end
toc=time()

#method3,usingreferencetothememberofcompositetypeobjects
#methodofchoice
tic=time()
for count in 1:1:iIterations
	for i=1:1:m
		this01=S01[i].M
		this02=S02[i].M
		this03=S03[i].M
		for j=1:1:n
			this03[j]=this01[j]+this02[j]
		end
	end
end
toc=time()