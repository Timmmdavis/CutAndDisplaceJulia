function DotLoopTest(x,y,xe,ye,a,Beta)

# Define material constant used in calculating influence coefficients
sb = sin(Beta); cb = cos(Beta);

R1= Array{Float64}(undef, length(x),1);
R2= Array{Float64}(undef, length(x),1);

for j = 1:length(x); 
	

    # Define array of local coordinates for the observation grid relative to
    #   the midpoint and orientation of the ith element
    # Refer to (Figure 56, C&S, p 91) and eqs 451 of C&S, p 57
    XB = (x[j]-xe)*cb + (y[j]-ye)*sb;
    YB = -(x[j]-xe)*sb + (y[j]-ye)*cb;

    # Calculate derivatives of the function f(x,y), eq 525 of C&S, p 81
    #   which are used to calculate the displacement and stress components
    # It is understood that X and Y refer to XB and YB
    # First abbreviate repeated terms in the derivatives of f(x,y):
    Y2 = YB^2;
    XMa = XB-a; XPa = XB+a;
    XMa2 = XMa^2; XPa2 = XPa^2;
    R1[j,1] = XMa2 + Y2; 
    R2[j,1] = XPa2 + Y2;
	
end
	
return(R1,R2)

end
