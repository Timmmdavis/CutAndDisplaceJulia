function AngDisDisp(x,y,z,cosA,sinA,bx,by,bz,E1,E2,cosA2,sinADE1)
# AngDisDisp calculates the "incomplete" displacements (without the 
# Burgers' function contribution) associated with an angular dislocation in
# an elastic full-space.

eta = y*cosA-z*sinA;
zeta = y*sinA+z*cosA;
r = sqrt(x^2+y^2+z^2);

# Avoid complex results for the logarithmic terms
if zeta>r
	zeta=r;
end
if z>r
	z = r;
end

b8p=bx/8/pi;
rMzeta=(r-zeta);
rMz=(r-z);

ux = b8p/E1*(x*y/r/rMz-x*eta/r/rMzeta);
vx = b8p/E1*(eta*sinA/rMzeta-y*eta/r/rMzeta+
    y.^2/r/rMz+E2*(cosA*log(rMzeta)-log(rMz)));
wx = b8p/E1*(eta*cosA/rMzeta-y/r-eta*z/r/rMzeta-
    E2*sinA*log(rMzeta));
uy = by/8/pi/E1*(x^2*cosA/r/rMzeta-x^2/r/rMz-
    E2*(cosA*log(rMzeta)-log(rMz)));
vy = by*x/8/pi/E1*(y*cosA/r/rMzeta-
    sinA*cosA/rMzeta-y/r/rMz);
wy = by*x/8/pi/E1*(z*cosA/r/rMzeta-
    cosA2/rMzeta+1/r);
	
uz = bz*sinADE1*(E2*log(rMzeta)-x.^2/r/rMzeta);
vz = bz*x*sinADE1*(sinA/rMzeta-y/r/rMzeta);
wz = bz*x*sinADE1*(cosA/rMzeta-z/r/rMzeta);

u = ux[1]+uy[1]+uz[1];
v = vx[1]+vy[1]+vz[1];
w = wx[1]+wy[1]+wz[1];

return(u,v,w)
end