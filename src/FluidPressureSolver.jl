###################################################################
# Fluid-filled cracks with a fixed volume (incompressible fluid)
###################################################################
#
# Called from SlipCalculator3D(...,BoundaryConditions::MixedBoundaryConditionsFluidVolume,...).
#
# Previously: guess a fluid pressure P, solve the friction LCP, measure the volume, and repeat
# (Brent/annealing) until the volume matched. That needs many LCP solves.
#
# Now: P is added as an extra unknown of the friction LCP, with one extra equation that fixes the
# volume. ONE solve gives the slip, the contact tractions AND the pressure.
#
# The friction LCP (Kaven et al. 2012, built in SlipCalculator3D(Scl,n,InvertedInfMatA,b,µ,Sf)) is:
#
#	[y] = [a][x] + [b]      with  x >= 0,  y >= 0,  x.*y = 0   (for each of the 5n rows)
#
#	where [y] = |+Dn |  and [x] = |-tn |      (contact tractions and slips, see SlipCalculator3D)
#	            |+Dss|            |-tss|
#	            |+Dds|            |-tds|
#	            |+Tss|            |-Dss|
#	            |+Tds|            |-Dds|
#
# A uniform fluid pressure P on crack k is the same as adding +P to Tn on its elements. This
# changes [b] linearly, [b] -> [b] + [g]*P. The column [g] holds the (scaled) Dn, Dss and Dds
# caused by a unit pressure on that crack. Its friction rows (4,5) are zero, as for the remote
# load in [b].
#
# The system we solve, one P and one volume equation per crack:
#
#	[y] = [a][x] + [b] + [g][P]       x >= 0, y >= 0, x.*y = 0      (5n complementarity rows)
#	[w]'[y]  = V                       sum over crack of Area*Dn       (1 volume row per crack)
#
#	where [w] has Area_i*Scl in the Dn rows of crack k (Dn = y[L1]*Scl) and zeros elsewhere.
#	P is free (no sign condition): it is whatever pressure holds the volume. So the tractions
#	given to SlipCalculator3D need no offset or shift.
#
# This is a mixed LCP. We solve it the same way as FischerNewton.fischer_newton: the
# complementarity rows are written as Fischer-Burmeister equations, and the volume rows are
# appended as plain equations:
#
#	phi(x_i,y_i) = sqrt(x_i^2+y_i^2) - x_i - y_i = 0      (zero exactly when x_i>=0, y_i>=0, x_i*y_i=0)
#	[w]'[y]/V - 1                                = 0      (volume, divided by V so it is O(1))
#
# Then Newton's method with a backtracking line search.
###################################################################


"""
	(g,w,CrackElements)=FluidVolumeSetup(C,Area,Scl,FractureFlag,n)

Build the extra column(s) [g] and volume row(s) [w] of the fluid volume problem (see the top of
FluidPressureSolver.jl). C is the inverted, scaled influence matrix (as in SlipCalculator3D).
One column of g and one column of w per crack.
"""
function FluidVolumeSetup(C,Area,Scl,FractureFlag,n)

	#Lengths for sub matricies/vectors (as SlipCalculator3D)
	L1=1:n; L2=n+1:2*n; L3=2*n+1:3*n

	#Elements belonging to each crack
	FractureFlag=vec(FractureFlag)
	NumOfFractures=convert(Int64,maximum(FractureFlag))
	CrackElements=[findall(FractureFlag.==k) for k=1:NumOfFractures]

	g=zeros(5*n,NumOfFractures) #Response of [y] to a unit pressure on each crack
	w=zeros(5*n,NumOfFractures) #Volume weights: [w]'[y] = volume of each crack
	for k=1:NumOfFractures
		Els=CrackElements[k]
		#Unit pressure = Tn of 1 on this crack's elements
		UnitTn=zeros(3*n)
		UnitTn[Els].=1.0
		#Scaled Dn, Dss, Dds that this causes (same as D=C*b in SlipCalculator3D).
		#Friction rows L4 and L5 stay 0.
		g[1:3*n,k]=C*UnitTn
		#Volume = sum(Area.*Dn) = sum(Area.*y[L1]*Scl)
		w[Els,k]=Area[Els].*Scl
	end

	return g,w,CrackElements
end


"""
	(x,y,P,Iterations,Converged)=FluidVolumeFischerNewton(a,b,g,w,V)

Solve the fluid volume mixed LCP (see the top of FluidPressureSolver.jl):
	y = a*x + b + g*P,   x >= 0, y >= 0, x.*y = 0,   w'*y = V
for x, y (5n) and the crack pressures P (one per crack).
"""
function FluidVolumeFischerNewton(a,b,g,w,V;Tolerance=1e-10,MaxIterations=100)

	N=length(b)            #5n complementarity rows
	m=length(V)            #one pressure / volume row per crack

	#Start: all elements open with zero contact traction (x=0), and the pressure that then gives
	#the right volume. With x=0, y=b+g*P, so w'*(b+g*P)=V.
	x=zeros(N)
	P=(w'*g)\(V-w'*b)

	#The complementarity residuals have units of the scaled displacements; divide them by a
	#typical size so they are comparable with the (dimensionless) volume residuals.
	Scale=max(maximum(abs.(b)),maximum(abs.(g*P)),1.0)

	#Residual F at (x,P): Fischer-Burmeister rows then volume rows
	function Residual(x,P)
		y=a*x+b+g*P
		Fcomp=(sqrt.(x.^2+y.^2)-x-y)./Scale
		Fvol=(w'*y)./V.-1.0
		return [Fcomp;Fvol],y
	end

	(F,y)=Residual(x,P)
	for Iterations=0:MaxIterations

		#Converged?
		if maximum(abs.(F))<=Tolerance
			return x,y,P,Iterations,true
		end
		if Iterations==MaxIterations
			break
		end

		#Jacobian of F wrt [x;P].
		#d(phi)/dx = x/sqrt(x^2+y^2)-1 = px,   d(phi)/dy = y/sqrt(x^2+y^2)-1 = qy,  and y=a*x+b+g*P
		#so  dFcomp/dx = (diag(px)+diag(qy)*a)/Scale,   dFcomp/dP = diag(qy)*g/Scale
		#    dFvol/dx  = (w'*a)./V,                      dFvol/dP  = (w'*g)./V
		#At x=y=0 phi has no derivative; px=qy=-1 is a valid choice there (as in fischer_newton).
		s=sqrt.(x.^2+y.^2)
		Zero=s.<1e-14*Scale
		px=ifelse.(Zero,-1.0,x./s.-1.0)
		qy=ifelse.(Zero,-1.0,y./s.-1.0)
		J=[(Diagonal(px)+Diagonal(qy)*a)./Scale   (Diagonal(qy)*g)./Scale;
		   (w'*a)./V                              (w'*g)./V              ]

		#Newton step (QR fallback: J can be singular when the tail first closes)
		Step=try
			J\(-F)
		catch err
			err isa SingularException || rethrow()
			qr(J,ColumnNorm())\(-F)
		end
		dx=Step[1:N]
		dP=Step[N+1:N+m]

		#Backtracking line search on 0.5*|F|^2, keeping x >= 0
		f=0.5*dot(F,F)
		τ=1.0
		while true
			xNew=max.(x+τ*dx,0.0)
			PNew=P+τ*dP
			(FNew,yNew)=Residual(xNew,PNew)
			if 0.5*dot(FNew,FNew)<=(1-2e-4*τ)*f || τ<1e-12
				x=xNew; P=PNew; F=FNew; y=yNew
				break
			end
			τ=τ/2
		end
	end

	return x,y,P,MaxIterations,false
end
