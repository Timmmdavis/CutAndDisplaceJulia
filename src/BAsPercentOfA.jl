function BAsPercentOfA(A,B)
# BAsPercentOfA - Compute B as a percentage of A. 

Percentage=(100.0./A).*B

#Puttings Infs to Nans (So we avoid this as max values)
Percentage[abs.(Percentage).==Inf].=NaN

return Percentage

end