function TensorTransformation3D(Txx1,Tyy1,Tzz1,Txy1,Txz1,Tyz1,A)
# TensorTransformation3D Transforms the coordinates of tensors,from x1y1z1 coordinate
# system to x2y2z2 coordinate system. "A" is the transformation matrix, 
# whose columns e1,e2 and e3 are the unit base vectors of the x1y1z1. The 
# coordinates of e1,e2 and e3 in A must be given in x2y2z2. The transpose 
# of A (i.e.., A') does the transformation from x2y2z2 into x1y1z1.

Txx2 = Array{Float64}(undef, length(Txx1),1); 
Tyy2 = Array{Float64}(undef, length(Txx1),1); 
Tzz2 = Array{Float64}(undef, length(Txx1),1); 
Txy2 = Array{Float64}(undef, length(Txx1),1); 
Txz2 = Array{Float64}(undef, length(Txx1),1); 
Tyz2 = Array{Float64}(undef, length(Txx1),1); 

for i=1:length(Txx1)
	Txx2[i] = A[1]^2*Txx1[i]+2*A[1]*A[4]*Txy1[i]+2*A[1]*A[7]*Txz1[i]+2*A[4]*A[7]*Tyz1[i]+A[4]^2*Tyy1[i]+A[7]^2*Tzz1[i];
	Tyy2[i] = A[2]^2*Txx1[i]+2*A[2]*A[5]*Txy1[i]+2*A[2]*A[8]*Txz1[i]+2*A[5]*A[8]*Tyz1[i]+A[5]^2*Tyy1[i]+A[8]^2*Tzz1[i];
	Tzz2[i] = A[3]^2*Txx1[i]+2*A[3]*A[6]*Txy1[i]+2*A[3]*A[9]*Txz1[i]+2*A[6]*A[9]*Tyz1[i]+A[6]^2*Tyy1[i]+A[9]^2*Tzz1[i];
	Txy2[i] = A[1]*A[2]*Txx1[i]+(A[1]*A[5]+A[2]*A[4])*Txy1[i]+(A[1]*A[8]+A[2]*A[7])*Txz1[i]+(A[8]*A[4]+A[7]*A[5])*Tyz1[i]+A[5]*A[4]*Tyy1[i]+A[7]*A[8]*Tzz1[i];
	Txz2[i] = A[1]*A[3]*Txx1[i]+(A[1]*A[6]+A[3]*A[4])*Txy1[i]+(A[1]*A[9]+A[3]*A[7])*Txz1[i]+(A[9]*A[4]+A[7]*A[6])*Tyz1[i]+A[6]*A[4]*Tyy1[i]+A[7]*A[9]*Tzz1[i];
	Tyz2[i] = A[2]*A[3]*Txx1[i]+(A[3]*A[5]+A[2]*A[6])*Txy1[i]+(A[3]*A[8]+A[2]*A[9])*Txz1[i]+(A[8]*A[6]+A[9]*A[5])*Tyz1[i]+A[5]*A[6]*Tyy1[i]+A[8]*A[9]*Tzz1[i];
end	
	
return(Txx2,Tyy2,Tzz2,Txy2,Txz2,Tyz2)

end