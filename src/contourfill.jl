function contourfill(x,y,v)
#A very basic version of contourf (like MATLAB)

#If you are getting errors do these before calling script{
##Get lengths (for reshapes later)
#dimx,dimy = size(x);

#x=reshape(x,dimx,dimy);
#y=reshape(y,dimx,dimy);
#v=reshape(v,dimx,dimy);
#}

#Get caxis limits
maxi=maximum(filter(!isnan, v));
mini=abs.(minimum(filter(!isnan, v)));
Top=maximum([maxi,mini])
#Steps from centre to top. 
steps=10; 
levels = [-Top:Top/steps:Top;]


#Close last fig
close() 
#draw contourf
contourf(x,y,v, levels=levels);

#add colourbar
cbar = colorbar()

end