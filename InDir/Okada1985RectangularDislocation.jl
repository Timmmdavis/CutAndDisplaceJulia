function Okada1985_RectangularDislocation(X,Y,Z,Strike,Dip,Length,Width,Rake,Ds,Dn,nu)

#OKADA85 Surface deformation due to a finite rectangular source.
#   #Tim comment: This function was originally called "OKADA85"
#	[uE,uN,uZ,uZE,uZN,uNN,uNE,uEN,uEE] = OKADA85(...
#	   E,N,DEPTH,STRIKE,DIP,LENGTH,WIDTH,RAKE,SLIP,OPEN)
#	computes displacements, tilts and strains at the surface of an elastic
#	half-space, due to a dislocation defined by RAKE, SLIP, and OPEN on a 
#	rectangular fault defined by orientation STRIKE and DIP, and size LENGTH and
#	WIDTH. The fault centroid is located (0,0,-DEPTH).
#
#	   E,N    : coordinates of observation points in a geographic referential 
#	            (East,North,Up) relative to fault centroid (units are described below)
#	   DEPTH  : depth of the fault centroid (DEPTH > 0)
#	   STRIKE : fault trace direction (0 to 360° relative to North), defined so 
#	            that the fault dips to the right side of the trace
#	   DIP    : angle between the fault and a horizontal plane (0 to 90°)
#	   LENGTH : fault length in the STRIKE direction (LENGTH > 0)
#	   WIDTH  : fault width in the DIP direction (WIDTH > 0)
#	   RAKE   : direction the hanging wall moves during rupture, measured relative
#	            to the fault STRIKE (-180 to 180°).
#	   SLIP   : dislocation in RAKE direction (length unit)
#	   OPEN   : dislocation in tensile component (same unit as SLIP)
#
#	returns the following variables (same matrix size as E and N):
#	   uN,uE,uZ        : displacements (unit of SLIP and OPEN)
#	   uZE,uZN         : tilts (in rad * FACTOR)
#	   uNN,uNE,uEN,uEE : horizontal strains POSITIVE = COMPRESSION (unit of FACTOR)
#
#	Length unit consistency: E, N, DEPTH, LENGTH, and WIDTH must have the same 
#	unit (e.g. km) which can be different from that of SLIP and OPEN (e.g. m) but
#	with a possible FACTOR on tilt and strain results (in this case, an 
#	amplification of km/m = 1000). To have FACTOR = 1 (tilt in radians and 
#	correct strain unit), use the same length unit for all aforesaid variables.


# Assigns input arguments
e = X;
n = Y;
depth = Z;
strike = Strike*pi/180;	# converting STRIKE in radian
dip = Dip*pi/180;	# converting DIP in radian ('delta' in Okada's equations)
L = Length;
W = Width;
rake = Rake*pi/180;	# converting RAKE in radian
slip = Ds;
U3 = Dn;
nnu;

# Defines dislocation in the fault plane system
U1 = cos(rake).*slip;
U2 = sin(rake).*slip;

# Converts fault coordinates (E,N,DEPTH) relative to centroid
# into Okada's reference system (X,Y,D)
d = depth + sin(dip).*W/2;	# d is fault's top edge
ec = e + cos(strike).*cos(dip).*W/2;
nc = n - sin(strike).*cos(dip).*W/2;
x = cos(strike).*nc + sin(strike).*ec + L/2;
y = sin(strike).*nc - cos(strike).*ec + cos(dip).*W;

# Variable substitution (independent from xi and eta)
p = y.*cos(dip) + d.*sin(dip);
q = y.*sin(dip) - d.*cos(dip);

#Displacements
ux = -U1/(2*pi) .* chinnery(ux_ss,x,p,L,W,q,dip,nu)  # strike-slip
	 -U2/(2*pi) .* chinnery(ux_ds,x,p,L,W,q,dip,nu)  # dip-slip
	 +U3/(2*pi) .* chinnery(ux_tf,x,p,L,W,q,dip,nu); # tensile fault

uy = -U1/(2*pi) .* chinnery(uy_ss,x,p,L,W,q,dip,nu)  # strike-slip
	 -U2/(2*pi) .* chinnery(uy_ds,x,p,L,W,q,dip,nu)  # dip-slip
	 +U3/(2*pi) .* chinnery(uy_tf,x,p,L,W,q,dip,nu); # tensile fault

uz = -U1/(2*pi) .* chinnery(uz_ss,x,p,L,W,q,dip,nu)  # strike-slip
	 -U2/(2*pi) .* chinnery(uz_ds,x,p,L,W,q,dip,nu)  # dip-slip
	 +U3/(2*pi) .* chinnery(uz_tf,x,p,L,W,q,dip,nu); # tensile fault
	
# Rotation from Okada's axes to geographic
ue = sin(strike).*ux - cos(strike).*uy;
un = cos(strike).*ux + sin(strike).*uy;
	
#Strain	
uxx = -U1/(2*pi) .* chinnery(uxx_ss,x,p,L,W,q,dip,nu)  # strike-slip
	  -U2/(2*pi) .* chinnery(uxx_ds,x,p,L,W,q,dip,nu)  # dip-slip
	  +U3/(2*pi) .* chinnery(uxx_tf,x,p,L,W,q,dip,nu); # tensile fault
uxy = -U1/(2*pi) .* chinnery(uxy_ss,x,p,L,W,q,dip,nu)  # strike-slip
	  -U2/(2*pi) .* chinnery(uxy_ds,x,p,L,W,q,dip,nu)  # dip-slip
	  +U3/(2*pi) .* chinnery(uxy_tf,x,p,L,W,q,dip,nu); # tensile fault
uyx = -U1/(2*pi) .* chinnery(uyx_ss,x,p,L,W,q,dip,nu)  # strike-slip
	  -U2/(2*pi) .* chinnery(uyx_ds,x,p,L,W,q,dip,nu)  # dip-slip
	  +U3/(2*pi) .* chinnery(uyx_tf,x,p,L,W,q,dip,nu); # tensile fault
uyy = -U1/(2*pi) .* chinnery(uyy_ss,x,p,L,W,q,dip,nu)  # strike-slip
	  -U2/(2*pi) .* chinnery(uyy_ds,x,p,L,W,q,dip,nu)  # dip-slip
	  +U3/(2*pi) .* chinnery(uyy_tf,x,p,L,W,q,dip,nu); # tensile fault
	  
# Rotation from Okada's axes to geographic
unn = cos(strike).^2*uxx + sin(2*strike).*(uxy + uyx)/2 + sin(strike).^2.*uyy;
une = sin(2*strike).*(uxx - uyy)/2 + sin(strike).^2.*uyx - cos(strike).^2.*uxy;
uen = sin(2*strike).*(uxx - uyy)/2 - cos(strike).^2.*uyx + sin(strike).^2.*uxy;
uee = sin(strike).^2*uxx - sin(2*strike).*(uyx + uxy)/2 + cos(strike).^2.*uyy;	

## Tilt
#uzx = -U1/(2*pi) .* chinnery(@uzx_ss,x,p,L,W,q,dip,nu)  # strike-slip
#	  -U2/(2*pi) .* chinnery(@uzx_ds,x,p,L,W,q,dip,nu)  # dip-slip
#	  +U3/(2*pi) .* chinnery(@uzx_tf,x,p,L,W,q,dip,nu); # tensile fault
#
#uzy = -U1/(2*pi) .* chinnery(@uzy_ss,x,p,L,W,q,dip,nu)  # strike-slip
#	  -U2/(2*pi) .* chinnery(@uzy_ds,x,p,L,W,q,dip,nu)  # dip-slip
#	  +U3/(2*pi) .* chinnery(@uzy_tf,x,p,L,W,q,dip,nu); # tensile fault
#
## Rotation from Okada's axes to geographic
#uze = -sin(strike).*uzx + cos(strike).*uzy;
#uzn = -cos(strike).*uzx - sin(strike).*uzy;pass functions into function julia

return(ux,uy,uz,unn,une,uen,uee)

end #end func	
	
	
##########################################################################

# Notes for I... and K... subfunctions:
#
#	1. original formulas use Lame's parameters as mu/(mu+lambda) which
#	   depends only on the Poisson's ratio = 1 - 2*nu
#	2. tests for cos(dip) == 0 are made with "cos(dip) > eps" 
#	   because cos(90*pi/180) is not zero but = 6.1232e-17


# =================================================================
# Chinnery's notation [equation (24) p. 1143]

# -----------------------------------------------------------------
function chinnery(f,x,p,L,W,q,dip,nu)
u = f(x,p,q,dip,nu)-f(x,p-W,q,dip,nu)- f(x-L,p,q,dip,nu)+ f(x-L,p-W,q,dip,nu);
return(u)
end


# =================================================================
# Displacement subfunctions

# strike-slip displacement subfunctions [equation (25) p. 1144]

# -----------------------------------------------------------------
function ux_ss(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = xi.*q./(R.*(R + eta)) ...
	+ I1(xi,eta,q,dip,nu,R).*sin(dip);
k = find(q~=0);
u(k) = u(k) + atan(xi(k).*eta(k)./(q(k).*R(k)));
return(u)
end

# -----------------------------------------------------------------
function uy_ss(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = (eta.*cos(dip) + q.*sin(dip)).*q./(R.*(R + eta)) ...
	+ q.*cos(dip)./(R + eta) ...
	+ I2(eta,q,dip,nu,R).*sin(dip);
return(u)
end

# -----------------------------------------------------------------
function uz_ss(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta.*sin(dip) - q.*cos(dip);
u = (eta.*sin(dip) - q.*cos(dip)).*q./(R.*(R + eta)) ...
	+ q.*sin(dip)./(R + eta) ...
	+ I4(db,eta,q,dip,nu,R).*sin(dip);
return(u)	
return(u)
end
	
# dip-slip displacement subfunctions [equation (26) p. 1144]
# -----------------------------------------------------------------
function ux_ds(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = q./R ...
	- I3(eta,q,dip,nu,R).*sin(dip).*cos(dip);
return(u)
end

# -----------------------------------------------------------------
function uy_ds(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = (eta.*cos(dip) + q.*sin(dip)).*q./(R.*(R + xi)) ...
	- I1(xi,eta,q,dip,nu,R).*sin(dip).*cos(dip);
k = find(q~=0);
u(k) = u(k) + cos(dip).*atan(xi(k).*eta(k)./(q(k).*R(k)));
return(u)
end

# -----------------------------------------------------------------
function uz_ds(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta.*sin(dip) - q.*cos(dip);
u = db.*q./(R.*(R + xi)) ...
	- I5(xi,eta,q,dip,nu,R,db).*sin(dip).*cos(dip);
k = find(q~=0);
u(k) = u(k) + sin(dip).*atan(xi(k).*eta(k)./(q(k).*R(k)));
return(u)
end

# tensile fault displacement subfunctions [equation (27) p. 1144]
# -----------------------------------------------------------------
function ux_tf(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = q.^2 ./(R.*(R + eta)) ...
	- I3(eta,q,dip,nu,R).*sin(dip).^2;
return(u)
end

# -----------------------------------------------------------------
function uy_tf(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = -(eta.*sin(dip) - q.*cos(dip)).*q./(R.*(R + xi)) ...
	- sin(dip).*xi.*q./(R.*(R + eta)) ...
	- I1(xi,eta,q,dip,nu,R).*sin(dip).^2;
k = find(q~=0);
u(k) = u(k) + sin(dip).*atan(xi(k).*eta(k)./(q(k).*R(k)));
return(u)
end

# -----------------------------------------------------------------
function uz_tf(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta.*sin(dip) - q.*cos(dip);
u = (eta.*cos(dip) + q.*sin(dip)).*q./(R.*(R + xi)) ...
	+ cos(dip).*xi.*q./(R.*(R + eta)) ...
	- I5(xi,eta,q,dip,nu,R,db).*sin(dip).^2;
k = find(q~=0);
u(k) = u(k) - cos(dip).*atan(xi(k).*eta(k)./(q(k).*R(k)));
return(u)
end

# I... displacement subfunctions [equations (28) (29) p. 1144-1145]
# -----------------------------------------------------------------
function I1(xi,eta,q,dip,nu,R)
db = eta.*sin(dip) - q.*cos(dip);
if cos(dip) > eps
	I = (1 - 2*nu) * (-xi./(cos(dip).*(R + db))) ...
		- sin(dip)./cos(dip).*I5(xi,eta,q,dip,nu,R,db);
else
	I = -(1 - 2*nu)/2 * xi.*q./(R + db).^2;
return(I)
end


# -----------------------------------------------------------------
function I2(eta,q,dip,nu,R)
I = (1 - 2*nu) * (-log(R + eta)) - I3(eta,q,dip,nu,R);
return(I)
end

# -----------------------------------------------------------------
function I3(eta,q,dip,nu,R)
yb = eta.*cos(dip) + q.*sin(dip);
db = eta.*sin(dip) - q.*cos(dip);
if cos(dip) > eps
	I = (1 - 2*nu) * (yb./(cos(dip)*(R + db)) - log(R + eta)) ...
		+ sin(dip)./cos(dip) * I4(db,eta,q,dip,nu,R);
else
	I = (1 - 2*nu)/2 * (eta./(R + db) + yb.*q./(R + db).^2 - log(R + eta));
return(I)
end

# -----------------------------------------------------------------
function I4(db,eta,q,dip,nu,R)
if cos(dip) > eps
	I = (1 - 2*nu) * 1./cos(dip) * (log(R + db) - sin(dip).*log(R + eta));
else
	I = -(1 - 2*nu) * q./(R + db);
return(I)
end

# -----------------------------------------------------------------
function I5(xi,eta,q,dip,nu,R,db)
X = sqrt(xi.^2 + q.^2);
if cos(dip) > eps
	I = (1 - 2*nu) * 2./cos(dip) ...
		.* atan((eta.*(X + q.*cos(dip)) + X.*(R + X).*sin(dip)) ...
			./(xi.*(R + X).*cos(dip)));
	I(xi==0) = 0;
else
	I = -(1 - 2*nu) * xi.*sin(dip)./(R + db);
return(I)
end

# =================================================================
# Tilt subfunctions

# strike-slip tilt subfunctions [equation (37) p. 1147]

# -----------------------------------------------------------------
function uzx_ss(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = -xi.*q.^2.*A(eta,R).*cos(dip) ...
	+ ((xi.*q)./R.^3 - K1(xi,eta,q,dip,nu,R)).*sin(dip);
return(u)
end

# -----------------------------------------------------------------
function uzy_ss(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta.*sin(dip) - q.*cos(dip);
yb = eta.*cos(dip) + q.*sin(dip);
u = (db.*q./R.^3).*cos(dip) ...
	+ (xi.^2.*q.*A(eta,R).*cos(dip) - sin(dip)./R + yb.*q./R.^3 ...
		- K2(xi,eta,q,dip,nu,R)).*sin(dip);
return(u)
end

# dip-slip tilt subfunctions [equation (38) p. 1147]

# -----------------------------------------------------------------
function uzx_ds(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta.*sin(dip) - q.*cos(dip);
u = db.*q./R.^3 ...
	+ q.*sin(dip)./(R.*(R + eta)) ...
	+ K3(xi,eta,q,dip,nu,R).*sin(dip).*cos(dip);
return(u)
end

# -----------------------------------------------------------------
function uzy_ds(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta.*sin(dip) - q.*cos(dip);
yb = eta.*cos(dip) + q.*sin(dip);
u = yb.*db.*q.*A(xi,R) ...
	- (2*db./(R.*(R + xi)) + xi.*sin(dip)./(R.*(R + eta))).*sin(dip) ...
	+ K1(xi,eta,q,dip,nu,R).*sin(dip).*cos(dip);
return(u)
end

# tensile fault tilt subfunctions [equation (39) p. 1147]

# -----------------------------------------------------------------
function uzx_tf(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = q.^2./R.^3.*sin(dip) ...
	- q.^3.*A(eta,R).*cos(dip) ...
	+ K3(xi,eta,q,dip,nu,R).*sin(dip).^2;
return(u)
end

# -----------------------------------------------------------------
function uzy_tf(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta.*sin(dip) - q.*cos(dip);
yb = eta.*cos(dip) + q.*sin(dip);
u = (yb.*sin(dip) + db.*cos(dip)).*q.^2.*A(xi,R) ...
	+ xi.*q.^2.*A(eta,R).*sin(dip).*cos(dip) ...
	- (2*q./(R.*(R + xi)) - K1(xi,eta,q,dip,nu,R)).*sin(dip).^2;
return(u)
end
	
# -----------------------------------------------------------------
function A(x,R)
a = (2*R + x)./(R.^3.*(R + x).^2);
return(a)
end

# K... tilt subfunctions [equations (40) (41) p. 1148]
# -----------------------------------------------------------------
function K1(xi,eta,q,dip,nu,R)
db = eta.*sin(dip) - q.*cos(dip);
if cos(dip) > eps
	K = (1 - 2*nu) * xi./cos(dip) .* (1./(R.*(R + db)) - sin(dip)./(R.*(R + eta)));
else
	K = (1 - 2*nu) * xi.*q./(R.*(R + db).^2);
return(K)
end


# -----------------------------------------------------------------
function K2(xi,eta,q,dip,nu,R)
K = (1 - 2*nu) * (-sin(dip)./R + q.*cos(dip)./(R.*(R + eta))) ...
	- K3(xi,eta,q,dip,nu,R);
return(K)
end

# -----------------------------------------------------------------
function K3(xi,eta,q,dip,nu,R)
db = eta.*sin(dip) - q.*cos(dip);
yb = eta.*cos(dip) + q.*sin(dip);
if cos(dip) > eps
	K = (1 - 2*nu) * 1./cos(dip) .* (q./(R.*(R + eta)) - yb./(R.*(R + db)));
else
	K = (1 - 2*nu) * sin(dip)./(R + db) .* (xi.^2./(R.*(R + db)) - 1);
return(K)
end


# =================================================================
# Strain subfunctions

# strike-slip strain subfunctions [equation (31) p. 1145]

# -----------------------------------------------------------------
function uxx_ss(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = xi.^2.*q.*A(eta,R) ...
	- J1(xi,eta,q,dip,nu,R).*sin(dip);
return(u)
end

# -----------------------------------------------------------------
function uxy_ss(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta.*sin(dip) - q.*cos(dip);
u = xi.^3.*db./(R.^3.*(eta.^2 + q.^2)) ...
	- (xi.^3.*A(eta,R) + J2(xi,eta,q,dip,nu,R)).*sin(dip);
return(u)
end
	
# -----------------------------------------------------------------
function uyx_ss(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = xi.*q./R.^3.*cos(dip) ...
	+ (xi.*q.^2.*A(eta,R) - J2(xi,eta,q,dip,nu,R)).*sin(dip);
return(u)
end

# -----------------------------------------------------------------
function uyy_ss(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
yb = eta.*cos(dip) + q.*sin(dip);
u = yb.*q./R.^3.*cos(dip) ...
	+ (q.^3.*A(eta,R).*sin(dip) - 2*q.*sin(dip)./(R.*(R + eta)) ...
		- (xi.^2 + eta.^2)./R.^3.*cos(dip) - J4(xi,eta,q,dip,nu,R)).*sin(dip);
return(u)
end	
	
# dip-slip strain subfunctions [equation (32) p. 1146]

# -----------------------------------------------------------------
function uxx_ds(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = xi.*q./R.^3 ...
	+ J3(xi,eta,q,dip,nu,R).*sin(dip).*cos(dip);
return(u)
end

# -----------------------------------------------------------------
function uxy_ds(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
yb = eta.*cos(dip) + q.*sin(dip);
u = yb.*q./R.^3 ...
	- sin(dip)./R ...
	+ J1(xi,eta,q,dip,nu,R).*sin(dip).*cos(dip);
return(u)
end

# -----------------------------------------------------------------
function uyx_ds(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
yb = eta.*cos(dip) + q.*sin(dip);
u = yb.*q./R.^3 ...
	+ q.*cos(dip)./(R.*(R + eta)) ...
	+ J1(xi,eta,q,dip,nu,R).*sin(dip).*cos(dip);
return(u)
end

# -----------------------------------------------------------------
function uyy_ds(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
yb = eta.*cos(dip) + q.*sin(dip);
u = yb.^2.*q.*A(xi,R) ...
	- (2*yb./(R.*(R + xi)) + xi.*cos(dip)./(R.*(R + eta))).*sin(dip) ...
	+ J2(xi,eta,q,dip,nu,R).*sin(dip).*cos(dip);
return(u)
end

# tensile fault strain subfunctions [equation (33) p. 1146]

# -----------------------------------------------------------------
function uxx_tf(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = xi.*q.^2.*A(eta,R) ...
	+ J3(xi,eta,q,dip,nu,R).*sin(dip).^2;
return(u)
end

# -----------------------------------------------------------------
function uxy_tf(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta.*sin(dip) - q.*cos(dip);
u = -db.*q./R.^3 ...
	- xi.^2.*q.*A(eta,R).*sin(dip) ...
	+ J1(xi,eta,q,dip,nu,R).*sin(dip).^2;
return(u)
end

# -----------------------------------------------------------------
function uyx_tf(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
u = q.^2./R.^3.*cos(dip) ...
	+ q.^3.*A(eta,R).*sin(dip) ...
	+ J1(xi,eta,q,dip,nu,R).*sin(dip).^2;
return(u)
end

# -----------------------------------------------------------------
function uyy_tf(xi,eta,q,dip,nu)
R = sqrt(xi.^2 + eta.^2 + q.^2);
db = eta.*sin(dip) - q.*cos(dip);
yb = eta.*cos(dip) + q.*sin(dip);
u = (yb.*cos(dip) - db.*sin(dip)).*q.^2.*A(xi,R) ...
	- q.*sin(2*dip)./(R.*(R + xi)) ...
	- (xi.*q.^2.*A(eta,R) - J2(xi,eta,q,dip,nu,R)).*sin(dip).^2;
return(u)
end

# J... tensile fault subfunctions [equations (34) (35) p. 1146-1147]
# -----------------------------------------------------------------
function J1(xi,eta,q,dip,nu,R)
db = eta.*sin(dip) - q.*cos(dip);
if cos(dip) > eps
	J = (1 - 2*nu) * 1./cos(dip) * (xi.^2./(R.*(R + db).^2) - 1./(R + db)) ...
		- sin(dip)./cos(dip)*K3(xi,eta,q,dip,nu,R);
else
	J = (1 - 2*nu)/2 * q./(R + db).^2 .* (2*xi.^2./(R.*(R + db)) - 1);
return(J)
end


# -----------------------------------------------------------------
function J2(xi,eta,q,dip,nu,R)
db = eta.*sin(dip) - q.*cos(dip);
yb = eta.*cos(dip) + q.*sin(dip);
if cos(dip) > eps
	J = (1 - 2*nu) * 1./cos(dip) * xi.*yb./(R.*(R + db).^2) ...
		- sin(dip)./cos(dip)*K1(xi,eta,q,dip,nu,R);
else
	J = (1 - 2*nu)/2 * xi.*sin(dip)./(R + db).^2 .* (2*q.^2./(R.*(R + db)) - 1);
return(J)
end

# -----------------------------------------------------------------
function J=J3(xi,eta,q,dip,nu,R)
J = (1 - 2*nu) * -xi./(R.*(R + eta)) ...
	- J2(xi,eta,q,dip,nu,R);
return(u)
end

# -----------------------------------------------------------------
function J4(xi,eta,q,dip,nu,R)
J = (1 - 2*nu) * (-cos(dip)./R - q.*sin(dip)./(R.*(R + eta))) ...
	- J1(xi,eta,q,dip,nu,R);
return(J)
end

