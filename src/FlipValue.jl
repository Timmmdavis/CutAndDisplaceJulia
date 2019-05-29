function FlipValue(a)
#Memoryless way of flipping a before putting into a func

for i=1:length(a)
	a[i]=-a[i]
end
return a
end