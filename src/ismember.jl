function ismember(a,b)
#a - large array
#b - value / row of values that we want to find in a 
#logical - returned logical list length of size(a,1)

#checking each col with loop
Total=zeros(size(a,1),1)
for i=1:size(a,2)
	Total.+=in.(a[:,i],[b[i]])
end
BinALogical=Total.==size(a,2)

return BinALogical
end