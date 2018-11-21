function AngDisStrain(x,y,z,cosA,sinA,bx,by,bz,nu,E1,E2,E3,cosA2,sinADE1)
# AngDisStrain calculates the strains associated with an angular 
# dislocation in an elastic full-space.

eta = y*cosA-z*sinA;
zeta = y*sinA+z*cosA;
E4=by*x/8/pi/E1;
E5=by/8/pi/E1;

x2 = x^2;
y2 = y^2;
z2 = z^2;
r2 = x2+y2+z2;
r = sqrt(r2);
r3 = r*r2;
rmz=(r-z);
rz = r*rmz;
r2z2 = r2*rmz^2;
r3z = r3*rmz;

W = zeta-r;
W2 = W^2;
Wr = W*r;
W2r = W2*r;
Wr3 = W*r3;
W2r2 = W2*r2;

C = (r*cosA-z)/Wr;
S = (r*sinA-y)/Wr;

# Partial derivatives of the Burgers' function
rFi_rx = (eta/r/(r-zeta)-y/r/rmz)/4/pi;
rFi_ry = (x/r/rmz-cosA*x/r/(r-zeta))/4/pi;
rFi_rz = (sinA*x/r/(r-zeta))/4/pi;

Exx = bx*(rFi_rx)+E3*(eta/Wr+eta*x2/W2r2-eta*x2/Wr3+y/rz-x2*y/r2z2-x2*y/r3z)-E4*((E2/Wr+x2/W2r2-x2/Wr3)*cosA+E2/rz-x2/r2z2-x2/r3z)+bz*x*sinADE1*(E2/Wr+x2/W2r2-x2/Wr3);

Eyy = by*(rFi_ry)+E3*((1/Wr+S^2-y2/Wr3)*eta+E2*y/rz-y^3/r2z2-y^3/r3z-2*nu*cosA*S)-E4*(1/rz-y2/r2z2-y2/r3z+(1/Wr+S^2-y2/Wr3)*cosA)+bz*x*sinADE1*(1/Wr+S^2-y2/Wr3);

Ezz = bz*(rFi_rz)+E3*(eta/W/r+eta*C^2-eta*z2/Wr3+y*z/r3+2*nu*sinA*C)-E4*((1/Wr+C^2-z2/Wr3)*cosA+z/r3)+bz*x*sinADE1*(1/Wr+C^2-z2/Wr3);

Exy = bx*(rFi_ry)/2+by*(rFi_rx)/2-E3*(x*y2/r2z2-nu*x/rz+x*y2/r3z-nu*x*cosA/Wr+eta*x*S/Wr+eta*x*y/Wr3)+E5*(x2*y/r2z2-nu*y/rz+x2*y/r3z+nu*cosA*S+x2*y*cosA/Wr3+x2*cosA*S/Wr)-bz*sinADE1*(nu*S+x2*S/Wr+x2*y/Wr3);

Exz = bx*(rFi_rz)/2+bz*(rFi_rx)/2-E3*(-x*y/r3+nu*x*sinA/Wr+eta*x*C/Wr+eta*x*z/Wr3)+E5*(-x2/r3+nu/r+nu*cosA*C+x2*z*cosA/Wr3+x2*cosA*C/Wr)-bz*sinADE1*(nu*C+x2*C/Wr+x2*z/Wr3);

Eyz = by*(rFi_rz)/2+bz*(rFi_ry)/2+E3*(y2/r3-nu/r-nu*cosA*C+nu*sinA*S+eta*sinA*cosA/W2-eta*(y*cosA+z*sinA)/W2r+eta*y*z/W2r2-eta*y*z/Wr3)-E4*(y/r3+sinA*cosA^2/W2-cosA*(y*cosA+z*sinA)/W2r+y*z*cosA/W2r2-y*z*cosA/Wr3)-bz*x*sinADE1*(y*z/Wr3-sinA*cosA/W2+(y*cosA+z*sinA)/W2r-y*z/W2r2);

return(Exx,Eyy,Ezz,Exy,Exz,Eyz)

end


