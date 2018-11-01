function meshgrid(x,y)
#Assuming vectors are both cols. Similar to MATLABS meshgrid func. 

dimx=length(x);
y=y';

#Like repmat
x=repeat(x, 1,length(y));
y=repeat(y,dimx,1);

return(x,y)

end
