function MatchingRow(A1,A2,i,j)
#Compare Els of two vectors
#memoryless version of: all(A1[i,:]==A2[j,:])
for p=1:6
    if A1[i,p]!=A2[j,p]
        return false
    end
end
return true
end

function MatchingRow(A1,A2,i,j,LengthofRows,flag)
#Compare Els of two vectors
#memoryless version of: all(A1[i,:]==A2[j,:])
for p=LengthofRows #assuming its assigned like LengthofRows=1:3
    if A1[i,p]!=A2[j,p]
    	flag=false #comes in as true
        return false
    end
end
return true
end