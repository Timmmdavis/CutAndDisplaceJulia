function meshgrid(x,y)
#Assuming vectors are both cols. Similar to MATLABS meshgrid func. 

dimy=length(y);
x=x';

#Like repmat
y=repeat(y, 1,length(x));
x=repeat(x,dimy,1);

return(x,y)

end
