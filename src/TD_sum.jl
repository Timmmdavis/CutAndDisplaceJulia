function TD_sum(ExxDn, EyyDn, EzzDn, ExyDn, ExzDn, EyzDn,
			    ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
				ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,
				UxDn,UyDn,UzDn,
				UxDss,UyDss,UzDss,
				UxDds,UyDds,UzDds)
#Simply sums the different modes of dislocations. (Safer than always writing this out)
#Could be written more elegantly 


if isempty(ExxDn)
	Exx=ExxDn;
	Eyy=EyyDn;
	Ezz=EzzDn;
	Exy=ExyDn;
	Exz=ExzDn;
	Eyz=EyzDn;
else
	#Accumulating arrays
	Exx=sum(ExxDn+ExxDss+ExxDds,dims=2);
	Eyy=sum(EyyDn+EyyDss+EyyDds,dims=2);
	Ezz=sum(EzzDn+EzzDss+EzzDds,dims=2);
	Exy=sum(ExyDn+ExyDss+ExyDds,dims=2);
	Exz=sum(ExzDn+ExzDss+ExzDds,dims=2);
	Eyz=sum(EyzDn+EyzDss+EyzDds,dims=2);
end

if isempty(UxDn)
	Ux=UxDn;
	Uz=UzDn;
	Uy=UyDn;
else
	Ux=sum(UxDn+UxDss+UxDds,dims=2);
	Uz=sum(UzDn+UzDss+UzDds,dims=2);
	Uy=sum(UyDn+UyDss+UyDds,dims=2);
end


return(Exx,Eyy,Ezz,Exy,Exz,Eyz,Ux,Uy,Uz)

end
