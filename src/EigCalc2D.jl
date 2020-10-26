function EigCalc3D(Sxx,Syy,Szz,Sxy,Sxz,Syz)
# EigCalc2d: Calculates the 2D principal stress/strain magnitudes and 
#                   directions from input tensors.
#   
# usage #1: For stress:
# [S1,S2,S1dir,S2dir]=EigCalc2d(Sxx,Syy,Sxy)
#
# usage #2: For strain:
# [E1,E2,E1dir,E2dir]=EigCalc2d(Exx,Eyy,Exy);
#
# Arguments: (input)
# Sxx,Syy,Sxy 		- The stress tensor components at a point, each can be a
#                    column vector.
#
# Arguments: (output)
# S1,S2       		- Principal stress component magnitudes Sigma 1 and
#                    Sigma 2.
#
# S1dir,S2dir       - Principal stress directions (direction cosines). Each
#                    will be a n*2 column vector [CosAx,CosAy] of this
#                    direction. 
#
# Example usage 1:
#
# #Calculating directions for a 2D stress tensor
# Sxx=0.2; Syy=-1.5; Sxy=1; 
# [S1,S2,S1dir,S2dir]=EigCalc2d(Sxx,Syy,Sxy)
# #Drawing
# DrawS1S2Directions( 0,0,S1dir,S2dir )
#
#  Author: Tim Davis
#  Copyright 2017, Tim Davis, Potsdam University\The University of Aberdeen



#Preallocating array
n = length(Sxx);
tensor=zeros(n*2,2);

#Filling this array with 3x3 tensors, every 3 rows is each point
#Indexing way to accumulate mats
N=2;
tensor[1:N:end,1] = Sxx[1:1:end,:];
tensor[1:N:end,2] = Sxy[1:1:end,:];
tensor[2:N:end,1] = Sxy[1:1:end,:];
tensor[2:N:end,2] = Syy[1:1:end,:];

#Eig can't handle nan's so we turn these to 0's and put the calculated s1s2s3 to nans after 
NanFlag = isnan.(tensor);
if any(NanFlag)
	tensor[NanFlag].=0;
end
#Calculating the eigen vectors and principal values (diagonal)
#Preallocating array
V=zeros(n*2,2); #Eigen vectors - http://eqseis.geosc.psu.edu/~cammon/HTML/UsingMATLAB/PDF/ML2#20Eigenvalues.pdf
D=zeros(n*2,2); #Eigen values -  also see Pollard 2005 addition resources MATLAB introduction WordDoc
lst=zeros(2)

for J=1:2:n*2

	lst,V[J:J+1,:]=LinearAlgebra.eigen(tensor[J:J+1,:])#,sortby= x -> -abs(x)
	D[J:J+1,:]=diagm(lst)
    V[J:J+1,:]=V[J:J+1,:]'; #Flipping the direction cosines as its easier to extract these like this. 

end


#Putting anywhere where there were nans in the tensor to nan
if any(NanFlag)
	D[NanFlag].=NaN;
	V[NanFlag].=NaN;
end
#Now getting a col vec of S2 S1 which is for each point
J=range(1,length(D[:,1]),step=2)
B=zeros(n*2)
for i = 1:n 
    B[J[i]:J[i]+1,:]=sort(diag(D[J[i]:J[i]+1,:]));
end  

#B output is S1S2S3 in a single column list, s1(a),s2(a),s3(a),s1(b),s2(b) etc. S1 S2 and S3 are then extracted to thier own arrays ready for 
#export from the function.
N=2;
S2 = B[1:N:end,:];
S1 = B[2:N:end,:];


#Also getting the direction cosines in the same manner. (note that X is first col, Y is second etc) 
S2dir = V[1:N:end,:];
S1dir = V[2:N:end,:];

return S1,S2,S1dir,S2dir
end