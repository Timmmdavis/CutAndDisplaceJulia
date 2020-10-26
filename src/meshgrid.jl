function meshgrid(x,y)
#Assuming vectors are both cols. Similar to MATLABS meshgrid func. 

dimx=length(x);
dimy=length(y);
x=reshape(x,1,dimx)

#Like repmat
x=repeat(x,dimy,1);
y=repeat(y,1,dimx);

return(x,y)

end


function meshgrid(x,y,z)
#In 3D. Similar to MATLABS meshgrid func. 

dimx=length(x);
dimy=length(y);
dimz=length(z);
x=reshape(x,1,dimx,1)
z=reshape(z,1,1,dimz)

#Like repmat
x=repeat(x, dimy,1,dimz);
y=repeat(y, 1,dimx,dimz);
z=repeat(z, dimy,dimx,1);
#z=reshape(z,dimx,dimy,dimz)

return(x,y,z)

end