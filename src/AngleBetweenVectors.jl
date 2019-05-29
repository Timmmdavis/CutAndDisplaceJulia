function AngleBetweenVectors(v1,v2,tmp,Ang)
#If calling inside loop preallocated tmp and Ang as 0. 
#outside of loop

#Check the input is a vector
if TestIsNotVec(v1)==true
	v1=vec(v1)
end
if TestIsNotVec(v2)==true
	v2=vec(v2)
end

#Assign to temporary array
tmp.=dot!(v1,v2)

#To avoid errors with acos
if abs(tmp)>1
	Ang.=0.
else 
	Ang.=acos(tmp)
end

return Ang
end

#memoryless vector check
function TestIsNotVec(a); 
	IsNotVec=typeof(a)!=Array{Int64,1} || typeof(a)!=Array{Float64,1}
end