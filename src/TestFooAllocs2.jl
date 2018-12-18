function TestFooAllocs2(Exx,Eyy,Ezz,Exy,Exz,Eyz,A,B)

# println([Exx[1],Eyy[1],Ezz[1],Exy[1],Exz[1],Eyz[1]]) #@info Exx[1] Eyy[1] Ezz[1] Exy[1] Exz[1] Eyz[1] 
# @time (exx,eyy,ezz,exy,exz,eyz)=TensorTransformation3D(Exx,Eyy,Ezz,Exy,Exz,Eyz,A)

# println([Exx[1],Eyy[1],Ezz[1],Exy[1],Exz[1],Eyz[1]])
# @time (Exx,Eyy,Ezz,Exy,Exz,Eyz)=TensorTransformation3D!(Exx,Eyy,Ezz,Exy,Exz,Eyz,A)



# if any(exx.!=Exx) || any(eyy.!=Eyy) || any(ezz.!=Ezz) || any(exy.!=Exy) || any(exz.!=Exz) || any(eyz.!=Eyz)
	# error("you have introduced an error")
# end


# @info Exx[1] Eyy[1] Ezz[1] Exy[1] Exz[1] Eyz[1] 
# @info exx[1] eyy[1] ezz[1] exy[1] exz[1] eyz[1] 



(Axx,Ayy,Azz,Axy,Axz,Ayz)=TensorTransformation3D(Exx,Eyy,Ezz,Exy,Exz,Eyz,A)
(exx,eyy,ezz,exy,exz,eyz)=TensorTransformation3D(Axx,Ayy,Azz,Axy,Axz,Ayz,B)


@info Exx[1] Eyy[1] Ezz[1] Exy[1] Exz[1] Eyz[1] 
@info exx[1] eyy[1] ezz[1] exy[1] exz[1] eyz[1] 


if any(!isequal(exx,Exx)) || any(!isequal(eyy,Eyy)) || any(!isequal(ezz,Ezz)) || any(!isequal(exy,Exy)) || any(!isequal(exz,Exz)) || any(!isequal(eyz,Eyz))
	error("you have introduced an error")
end



end