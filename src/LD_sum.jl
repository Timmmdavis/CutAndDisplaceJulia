function LD_sum(SxxDs,SxxDn,SyyDs,
				SyyDn,SxyDs,SxyDn,
				UxDs,UxDn,
				UyDs,UyDn)
#Simply sums the different modes of dislocations. (Safer than always writing this out)
#Could be written more elegantly 

#Accumulating arrays
Sxx=sum(SxxDs.+SxxDn,dims=2);
Syy=sum(SyyDs.+SyyDn,dims=2);
Sxy=sum(SxyDs.+SxyDn,dims=2);
Ux=sum(UxDs.+UxDn,dims=2);
Uy=sum(UyDs.+UyDn,dims=2);

return(Sxx,Syy,Sxy,Ux,Uy)

end
