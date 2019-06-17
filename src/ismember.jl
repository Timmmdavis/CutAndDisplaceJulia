function ismember(a,b)
#a - large array
#b - value / row of values that we want to find in a 
#logical - returned logical list length of size(a,1)

#checking each col with loop
Total=zeros(size(a,1))
n=size(a,1)
n2=size(a,2)
for j=1:n
	for i=1:n2
		Total[j]+=in(a[j,i],b[i])
	end
end
BinALogical=vec(Total.==n2)

return BinALogical
end

function ismember(a,b,n::Int,n2::Int)
#a - large array
#b - value / row of values that we want to find in a 
#logical - returned logical list length of size(a,1)
#n - size(a,1)
#n - size(a,2)

#checking each col with loop
Total=fill(0,size(a,1))
Flag=fill(false,size(a,1))
for j=1:n
	for i=1:n2
		Total[j]+=in(a[j,i],b[i])
	end
	#convert to Int
	if Total[j]==n2
		Flag[j]=true
	end
end

return Flag
end