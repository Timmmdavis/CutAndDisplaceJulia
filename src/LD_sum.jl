function LD_sum(SxxDs,SxxDn,SyyDs,SyyDn,SxyDs,SxyDn,UxDs,UxDn,UyDs,UyDn)
#Simply sums the different modes of dislocations. (Safer than always writing this out)

#Accumulating arrays
Sxx=SxxDs.+SxxDn;
Syy=SyyDs.+SyyDn;
Sxy=SxyDs.+SxyDn;
Ux=UxDs.+UxDn;
Uy=UyDs.+UyDn;

return(Sxx,Syy,Sxy,Ux,Uy)

end
