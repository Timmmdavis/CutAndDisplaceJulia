
# TensorTransformation3D Transforms the coordinates of tensors,from x1y1z1 coordinate
# system to x2y2z2 coordinate system "A" is the transformation matrix, 
# whose columns e1,e2 and e3 are the unit base vectors of the x1y1z The 
# coordinates of e1,e2 and e3 in A must be given in x2y2z2 The transpose 
# of A (ie, A') does the transformation from x2y2z2 into x1y1z


@inline function TensorTransformation3D!(Txx::Array,Tyy::Array,Tzz::Array,Txy::Array,Txz::Array,Tyz::Array,A)
@fastmath @simd for i=1:length(Txx)

	txx = A[1]^2 *Txx[i] +2*A[1]*A[4]*Txy[i] +2*A[1]*A[7]*Txz[i] +2*A[4]*A[7]*Tyz[i]+A[4]^2 *Tyy[i]+A[7]^2 *Tzz[i];
	
	tyy = A[2]^2 *Txx[i] +2*A[2]*A[5]*Txy[i] +2*A[2]*A[8]*Txz[i] +2*A[5]*A[8]*Tyz[i]+A[5]^2 *Tyy[i]+A[8]^2 *Tzz[i];
	
	tzz = A[3]^2 *Txx[i] +2*A[3]*A[6]*Txy[i] +2*A[3]*A[9]*Txz[i] +2*A[6]*A[9]*Tyz[i]+A[6]^2 *Tyy[i]+A[9]^2 *Tzz[i];
	
	txy = A[1]*A[2]*Txx[i]+(A[1]*A[5]+A[2]*A[4])*Txy[i]+(A[1]*A[8]+A[2]*A[7])*Txz[i]+(A[4]*A[8]+A[5]*A[7])*Tyz[i]+A[4]*A[5]*Tyy[i]+A[7]*A[8]*Tzz[i];
	
	txz = A[1]*A[3]*Txx[i]+(A[1]*A[6]+A[3]*A[4])*Txy[i]+(A[1]*A[9]+A[3]*A[7])*Txz[i]+(A[4]*A[9]+A[6]*A[7])*Tyz[i]+A[4]*A[6]*Tyy[i]+A[7]*A[9]*Tzz[i];
	
	tyz = A[2]*A[3]*Txx[i]+(A[3]*A[5]+A[2]*A[6])*Txy[i]+(A[3]*A[8]+A[2]*A[9])*Txz[i]+(A[6]*A[8]+A[5]*A[9])*Tyz[i]+A[5]*A[6]*Tyy[i]+A[8]*A[9]*Tzz[i];
	
	Txx[i]=txx;
	Tyy[i]=tyy;
	Tzz[i]=tzz;
	Txy[i]=txy;
	Txz[i]=txz;
	Tyz[i]=tyz;
end
return(Txx,Tyy,Tzz,Txy,Txz,Tyz)
end

@inline function TensorTransformation3D!(Txx::Float64,Tyy::Float64,Tzz::Float64,Txy::Float64,Txz::Float64,Tyz::Float64,A)


	txx = A[1]^2 *Txx +2*A[1]*A[4]*Txy +2*A[1]*A[7]*Txz +2*A[4]*A[7]*Tyz+A[4]^2 *Tyy+A[7]^2 *Tzz;
	tyy = A[2]^2 *Txx +2*A[2]*A[5]*Txy +2*A[2]*A[8]*Txz +2*A[5]*A[8]*Tyz+A[5]^2 *Tyy+A[8]^2 *Tzz;
	tzz = A[3]^2 *Txx +2*A[3]*A[6]*Txy +2*A[3]*A[9]*Txz +2*A[6]*A[9]*Tyz+A[6]^2 *Tyy+A[9]^2 *Tzz;
	txy = A[1]*A[2]*Txx+(A[1]*A[5]+A[2]*A[4])*Txy+(A[1]*A[8]+A[2]*A[7])*Txz+(A[4]*A[8]+A[5]*A[7])*Tyz+A[4]*A[5]*Tyy+A[7]*A[8]*Tzz;
	txz = A[1]*A[3]*Txx+(A[1]*A[6]+A[3]*A[4])*Txy+(A[1]*A[9]+A[3]*A[7])*Txz+(A[4]*A[9]+A[6]*A[7])*Tyz+A[4]*A[6]*Tyy+A[7]*A[9]*Tzz;
	tyz = A[2]*A[3]*Txx+(A[3]*A[5]+A[2]*A[6])*Txy+(A[3]*A[8]+A[2]*A[9])*Txz+(A[6]*A[8]+A[5]*A[9])*Tyz+A[5]*A[6]*Tyy+A[8]*A[9]*Tzz;
	
	Txx=txx;
	Tyy=tyy;
	Tzz=tzz;
	Txy=txy;
	Txz=txz;
	Tyz=tyz;

return(Txx,Tyy,Tzz,Txy,Txz,Tyz)
end