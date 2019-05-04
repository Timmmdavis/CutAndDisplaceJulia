function ismember(a,b)
#a - large array
#b - value / row of values that we want to find in a 
#logical - returned logical list length of size(a,1)

BinALogical=sum(Int.(in.(a,[b])),dims=2).==size(a,2)

return BinALogical
end