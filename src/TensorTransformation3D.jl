function TensorTransformation3D(Txx1,Tyy1,Tzz1,Txy1,Txz1,Tyz1,A)
# TensorTransformation3D Transforms the coordinates of tensors,from x1y1z1 coordinate
# system to x2y2z2 coordinate system. "A" is the transformation matrix, 
# whose columns e1,e2 and e3 are the unit base vectors of the x1y1z1. The 
# coordinates of e1,e2 and e3 in A must be given in x2y2z2. The transpose 
# of A (i.e.., A') does the transformation from x2y2z2 into x1y1z1.

Txx2 = A[1].^2 .*Txx1.+2 .*A[1].*A[4].*Txy1.+2 .*A[1].*A[7].*Txz1.+2 .*A[4].*A[7].*Tyz1.+A[4].^2 .*Tyy1.+A[7].^2 .*Tzz1;
Tyy2 = A[2].^2 .*Txx1.+2 .*A[2].*A[5].*Txy1.+2 .*A[2].*A[8].*Txz1.+2 .*A[5].*A[8].*Tyz1.+A[5].^2 .*Tyy1.+A[8].^2 .*Tzz1;
Tzz2 = A[3].^2 .*Txx1.+2 .*A[3].*A[6].*Txy1.+2 .*A[3].*A[9].*Txz1.+2 .*A[6].*A[9].*Tyz1.+A[6].^2 .*Tyy1.+A[9].^2 .*Tzz1;
Txy2 = A[1].*A[2].*Txx1.+(A[1].*A[5].+A[2].*A[4]).*Txy1.+(A[1].*A[8].+A[2].*A[7]).*Txz1.+(A[8].*A[4].+A[7].*A[5]).*Tyz1.+A[5].*A[4].*Tyy1.+A[7].*A[8].*Tzz1;
Txz2 = A[1].*A[3].*Txx1.+(A[1].*A[6].+A[3].*A[4]).*Txy1.+(A[1].*A[9].+A[3].*A[7]).*Txz1.+(A[9].*A[4].+A[7].*A[6]).*Tyz1.+A[6].*A[4].*Tyy1.+A[7].*A[9].*Tzz1;
Tyz2 = A[2].*A[3].*Txx1.+(A[3].*A[5].+A[2].*A[6]).*Txy1.+(A[3].*A[8].+A[2].*A[9]).*Txz1.+(A[8].*A[6].+A[9].*A[5]).*Tyz1.+A[5].*A[6].*Tyy1.+A[8].*A[9].*Tzz1;


# Txx2=zeros(size(Txx1));
# Tyy2=zeros(size(Txx1));
# Tzz2=zeros(size(Txx1));
# Txy2=zeros(size(Txx1));
# Txz2=zeros(size(Txx1));
# Tyz2=zeros(size(Txx1));
# for i=1:length(Txx2)
    # Quat=A;
    # Tensor=[[Txx1[i,:] Txy1[i,:] Txz1[i,:]]; [Txy1[i,:] Tyy1[i,:] Tyz1[i,:]]; [Txz1[i,:] Tyz1[i,:] Tzz1[i,:]]];
    # CartStress=Quat*Tensor*Quat';
    # Txx2[i]=CartStress[1,1];
    # Tyy2[i]=CartStress[2,2];
    # Tzz2[i]=CartStress[3,3];
    # Txy2[i]=CartStress[1,2];
    # Txz2[i]=CartStress[1,3];
    # Tyz2[i]=CartStress[2,3];
# end


return(Txx2,Tyy2,Tzz2,Txy2,Txz2,Tyz2)

end