function cart2pol(x,y)
# cart2pol: hash of MATLABS cart2pol function

#2D mapping
θ=atan.(y,x);
ρ=sqrt.(x.^2+y.^2);

return(θ,ρ) 

end

function cart2pol(x,y,z)
# cart2pol: hash of MATLABS cart2pol function

#3D mapping
(θ,ρ)=cart2pol(x,y)
z=z;

return(θ,ρ,z) 

end