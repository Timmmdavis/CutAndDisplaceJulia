function pol2cart(θ,ρ)
# pol2cart: hash of MATLABS pol2cart function

#2D mapping
x=ρ.*cos.(θ)
y=ρ.*sin.(θ)

return(x,y) 

end

function pol2cart(θ,ρ,z)
# pol2cart: hash of MATLABS pol2cart function

#3D mapping
(x,y)=pol2cart(θ,ρ)
z=z;

return(x,y,z) 

end