###################################################################
# Fluid-filled cracks with a fixed volume (incompressible fluid)
###################################################################
#
# Called from SlipCalculator3D(...,BoundaryConditions::MixedBoundaryConditionsFluidVolume,...).
#
# The fluid pressure P of each crack is an extra unknown of the friction LCP, with one extra equation
# per crack that fixes its volume. ONE solve gives the slip, the contact tractions AND the pressure.
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
#	P is free (no sign condition): it is whatever pressure holds the volume.
#
# The complementarity rows are written as Fischer-Burmeister equations and the volume rows are
# appended as plain equations, then solved with Newton's method and a backtracking line search:
#
#	phi(x_i,y_i) = sqrt(x_i^2+y_i^2) - x_i - y_i = 0      (zero exactly when x_i>=0, y_i>=0, x_i*y_i=0)
#	[w]'[y]/V - 1                                = 0      (volume, divided by V so it is O(1))
#
# Speed: at an open element (y>0, x=0) the Fischer-Burmeister row only depends on x_i, so its Newton step
# is explicit. Only the remaining rows (contact / sliding) and the volume rows need a linear solve: the
# "reduced" system, typically half the size, so each LU is several times cheaper.
#
# Singular Jacobians: on an open element the friction strength mu*0+Sf is zero, so the upper and lower
# friction bounds coincide; only the difference y2-x4 (=Dss) is fixed and x4, y2 can move together.
# The Newton Jacobian is therefore usually exactly singular. This is a property of the formulation, not
# bad scaling. It is handled by a tiny diagonal shift (cheap), and where the linearised system is also
# inconsistent (the shifted step blows up and the line search fails) by a pivoted QR, which gives the
# minimum-norm least-squares step. If the reduced solve still does not converge, the full system is
# solved with pivoted QR as a fallback (slow but robust).
###################################################################


"""
	(g,w,CrackElements)=FluidVolumeSetup(C,Area,Scl,FractureFlag,n)

Build the extra column(s) [g] and volume row(s) [w] of the fluid volume problem (see the top of
FluidPressureSolver.jl). C is the inverted, scaled influence matrix (as in SlipCalculator3D).
One column of g and one column of w per crack.
"""
function FluidVolumeSetup(C,Area,Scl,FractureFlag,n)

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
	(x,y,P,Iterations,Converged)=FluidVolumeFischerNewton(a,b,g,w,V;Tolerance=1e-10,MaxIterations=100)

Solve the fluid volume mixed LCP (see the top of FluidPressureSolver.jl):
	y = a*x + b + g*P,   x >= 0, y >= 0, x.*y = 0,   w'*y = V
for x, y (5n) and the crack pressures P (one per crack).

Uses the reduced Newton system. It normally converges in 7-15 iterations, so it is stopped after 30, or
earlier if it stalls, and the full Newton system (pivoted QR) is then solved instead.
"""
function FluidVolumeFischerNewton(a,b,g,w,V;Tolerance=1e-10,MaxIterations=100)
	out=FluidVolumeNewton(a,b,g,w,V,true;Tolerance=Tolerance,MaxIterations=min(MaxIterations,30))
	out[5] && return out
	@warn "Fluid volume solve: reduced Newton did not converge ($(out[4]) iterations); solving the full system"
	return FluidVolumeNewton(a,b,g,w,V,false;Tolerance=Tolerance,MaxIterations=MaxIterations)
end

"""
	(x,y,P,Iterations,Converged)=FluidVolumeNewton(a,b,g,w,V,Reduced;Tolerance=1e-10,MaxIterations=100)

Newton's method for the fluid volume mixed LCP (used by FluidVolumeFischerNewton). Reduced=true solves the
reduced system and stops early if it stalls; Reduced=false solves the full system.
"""
function FluidVolumeNewton(a,b,g,w,V,Reduced;Tolerance=1e-10,MaxIterations=100)

	N=length(b)            #5n complementarity rows
	m=length(V)            #one pressure / volume row per crack

	#Start: all elements open with zero contact traction (x=0), and the pressure that then gives
	#the right volume. With x=0, y=b+g*P, so w'*(b+g*P)=V.
	x=zeros(N)
	P=(w'*g)\(V-w'*b)

	#The complementarity residuals have units of the scaled displacements; divide them by a
	#typical size so they are comparable with the (dimensionless) volume residuals.
	Scale=max(maximum(abs,b),maximum(abs,g*P),1.0)

	#Volume rows of the Jacobian (constant)
	wa=(w'*a)./V
	wg=(w'*g)./V

	#Residual F at (x,P): Fischer-Burmeister rows then volume rows
	function Residual(x,P)
		y=a*x+b+g*P
		Fcomp=(sqrt.(x.^2+y.^2)-x-y)./Scale
		Fvol=(w'*y)./V.-1.0
		return [Fcomp;Fvol],y
	end

	Singular=false         #true once an LU has found the Jacobian singular (it then stays singular)
	History=Float64[]      #residual per iteration, to spot a stalled solve
	(F,y)=Residual(x,P)
	for Iterations=0:MaxIterations

		normF=maximum(abs,F)
		normF<=Tolerance && return x,y,P,Iterations,true
		Iterations==MaxIterations && break
		push!(History,normF)
		#Reduced solve stalled (residual not 10x smaller than 8 iterations ago): give up, the caller falls back
		Reduced && length(History)>8 && normF>History[end-8]/10 && return x,y,P,Iterations,false

		#Jacobian of F wrt [x;P].
		#d(phi)/dx = x/sqrt(x^2+y^2)-1 = px,   d(phi)/dy = y/sqrt(x^2+y^2)-1 = qy,  and y=a*x+b+g*P
		#so  dFcomp/dx = (diag(px)+diag(qy)*a)/Scale,   dFcomp/dP = diag(qy)*g/Scale
		#    dFvol/dx  = (w'*a)./V,                      dFvol/dP  = (w'*g)./V
		#At x=y=0 phi has no derivative; px=qy=-1 is a valid choice there (as in fischer_newton).
		s=sqrt.(x.^2+y.^2)
		Zero=s.<1e-14*Scale
		px=ifelse.(Zero,-1.0,x./s.-1.0)
		qy=ifelse.(Zero,-1.0,y./s.-1.0)

		if Reduced
			#Open elements (qy=0): Jacobian row is px_i/Scale on the diagonal only, so dx_i is explicit.
			Triv=findall(abs.(qy).<1e-9)
			Act=findall(abs.(qy).>=1e-9)
			na=length(Act)
			dxTriv=-F[Triv].*Scale./px[Triv]
			J=Matrix{Float64}(undef,na+m,na+m)
			@views J[1:na,1:na].=qy[Act].*a[Act,Act]./Scale
			for (k,i) in enumerate(Act)
				J[k,k]+=px[i]/Scale
			end
			@views J[1:na,na+1:end].=qy[Act].*g[Act,:]./Scale
			@views J[na+1:end,1:na].=wa[:,Act]
			J[na+1:end,na+1:end].=wg
			rhs=[-F[Act]-qy[Act].*(a[Act,Triv]*dxTriv)./Scale; -F[N+1:N+m]-wa[:,Triv]*dxTriv]
			(sol,Singular,Regularised)=NewtonSolve(J,rhs,Singular,true)
			Expand=sol->(dx=zeros(N); dx[Triv]=dxTriv; dx[Act]=sol[1:na]; (dx,sol[na+1:na+m]))
		else
			J=[(Diagonal(px)+Diagonal(qy)*a)./Scale   (Diagonal(qy)*g)./Scale;
			   wa                                       wg                    ]
			(sol,Singular,Regularised)=NewtonSolve(J,-F,Singular,false)
			Expand=sol->(sol[1:N],sol[N+1:N+m])
		end
		(dx,dP)=Expand(sol)

		#Backtracking line search on 0.5*|F|^2, keeping x >= 0.
		#If a regularised step needs a tiny τ, the linearised system was inconsistent (not just singular) and the
		#diagonal shift has blown the step up: redo this step with pivoted QR (minimum-norm least squares).
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
			if Regularised && τ<2.0^-10
				(dx,dP)=Expand(qr(J,ColumnNorm())\(Reduced ? rhs : -F))
				Regularised=false
				τ=1.0
			end
		end
	end

	return x,y,P,MaxIterations,false
end

"""
	(s,Singular,Regularised)=NewtonSolve(J,r,Singular,Regularise)

Newton step s from J*s=r. A plain LU is tried until the Jacobian has once been found singular (it then stays
singular, see the top of FluidPressureSolver.jl). For a singular J: Regularise=true uses an LU of J plus a tiny
diagonal shift (1e-8 of its largest entry: picks the minimum-change step in the undetermined directions; cheap);
otherwise, or if that also fails, a pivoted QR (minimum-norm least squares; robust, ~10x slower than LU).
Regularised says whether the shifted LU was used. J is not changed.
"""
function NewtonSolve(J,r,Singular,Regularise)
	if !Singular
		F=lu(J,check=false)
		issuccess(F) && return F\r,false,false
		Singular=true
	end
	if Regularise
		F=lu(J+(1e-8*maximum(abs,J))*I,check=false)
		issuccess(F) && return F\r,true,true
	end
	return qr(J,ColumnNorm())\r,true,false
end
