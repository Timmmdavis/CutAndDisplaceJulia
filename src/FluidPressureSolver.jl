# Fluid-filled cracks solved in ONE friction-LCP solve: the fluid pressure of each crack is an extra
# unknown and its volume (or mass) an extra equation. Used by
#   SlipCalculator3D(..., ::MixedBoundaryConditionsFluidVolume, ...)   incompressible, fixed volume
#   SlipCalculator3D(..., ::MixedBoundaryConditionsFluidMass, ...)     compressible, fixed mass
#
# Friction LCP: Kaven et al. (2012), as in SlipCalculator3D(Scl,n,InvertedInfMatA,b,µ,Sf):
#   y = a*x + b0 + G*u,  x >= 0, y >= 0, x.*y = 0   (a, b0: 5n friction LCP; u: fluid unknowns)
# Fluid unknowns per crack k: P_k (pressure at elevation zref_k) and s_k = ρ̄_k/ρ0_k.
#   fluid pressure on element i of crack k:  p_i = P_k + s_k*ρ0_k*g*(zref_k - z_i)
#   equation of state (linearised):          ρ = ρ0*(1 + (p_i + pshift - pa)/K)
# Extra equations per crack:
#   r1: Σ A_i Dn_i (1 + (p_i + pshift - pa)/K) = V0      (mass conserved; V0 = mass/ρ0)
#   r2: s * Σ A_i Dn_i                          = V0      (ρ̄ * volume = mass)
# Incompressible (fixed volume): ρ0 = 0 (fluid weight, if any, is in the tractions you supply), K = Inf.
# Isothermal ideal gas: K = pa = reference pressure, ρ0 = density at pa (the linear EOS is then exact).
#
# Solved with Fischer-Burmeister semismooth Newton (as FischerNewton.fischer_newton) and a projected
# Armijo line search; dense direct solves for the Newton step.

struct FluidCrack
    idx::Vector{Int}   # element indices of this crack
    V0::Float64        # volume to hold (incompressible) or mass/ρ0 (compressible)
    ρ0::Float64        # fluid density at pa (0 for incompressible: no internal weight)
    K::Float64         # fluid bulk modulus (Inf = incompressible)
    zref::Float64      # elevation where P is defined (top midpoint of the crack)
end

function MakeFluidCracks(FractureFlag, MidPoint, V0, ρ0, K)
    Flag = vec(FractureFlag)
    nf = convert(Int, maximum(Flag))
    pick(v, k) = v isa Number ? Float64(v) : Float64(v[k])
    cracks = FluidCrack[]
    for k in 1:nf
        idx = findall(Flag .== k)
        push!(cracks, FluidCrack(idx, pick(V0, k), pick(ρ0, k), pick(K, k), maximum(MidPoint[idx, 3])))
    end
    return cracks
end

# Scaled displacement response (3n × 2nf) to unit pressure / unit weight column on each crack
function FluidColumns(Ainv, MidPoint, cracks, n; g = 9.81)
    nf = length(cracks)
    z = MidPoint[:, 3]
    Gm = zeros(3n, 2nf)
    AinvTn = view(Ainv, :, 1:n)
    for (k, c) in enumerate(cracks)
        unitP = zeros(n); unitP[c.idx] .= 1.0
        grad = zeros(n); grad[c.idx] .= c.ρ0 * g .* (c.zref .- z[c.idx])
        Gm[:, k] = AinvTn * unitP
        Gm[:, nf + k] = AinvTn * grad
    end
    return Gm
end

# Mass/volume equations r(y,u), with derivatives wrt y (Dn rows only) and u
function FluidRows(cracks, Area, MidPoint, Scl, N; g = 9.81, pshift = 0.0, pa = 0.0)
    nf = length(cracks)
    z = MidPoint[:, 3]
    A = vec(Area)
    function rows(y, u)
        r = zeros(2nf); Ry = zeros(2nf, N); Ru = zeros(2nf, 2nf)
        for (k, c) in enumerate(cracks)
            P = u[k]; s = u[nf + k]
            V = 0.0; Vp = 0.0; Vh = 0.0
            for i in c.idx
                Dn = Scl * y[i]; dV = A[i] * Dn
                h = c.ρ0 * g * (c.zref - z[i])
                p = P + s * h + pshift - pa
                V += dV; Vp += dV * (1 + p / c.K); Vh += dV * h
                Ry[k, i] = Scl * A[i] * (1 + p / c.K) / c.V0
                Ry[nf + k, i] = s * Scl * A[i] / c.V0
            end
            r[k] = Vp / c.V0 - 1
            r[nf + k] = s * V / c.V0 - 1
            Ru[k, k] = V / c.K / c.V0
            Ru[k, nf + k] = Vh / c.K / c.V0
            Ru[nf + k, nf + k] = V / c.V0
        end
        return r, Ry, Ru
    end
    return rows
end

# Starting guess: incompressible, no-contact pressures that give V = V0 (s = 1)
function FluidLinearGuess(bn, Gn, cracks, Area, Scl)
    nf = length(cracks); A = vec(Area)
    Wt = zeros(nf, length(bn))
    for (k, c) in enumerate(cracks)
        Wt[k, c.idx] .= Scl .* A[c.idx]
    end
    rhs = [c.V0 for c in cracks] .- Wt * (bn .+ Gn[:, nf+1:2nf] * ones(nf))
    P = (Wt * Gn[:, 1:nf]) \ rhs
    return vcat(P, ones(nf))
end

"""
    (x, u, y, iters, ok) = BorderedFB(a, b0, G5, rows, u0; x0, rtol, maxit)

Friction LCP y = a*x + b0 + G5*u plus fluid rows r(y,u) = 0, by Fischer-Burmeister Newton.
Residuals of the complementarity rows are divided by a traction scale so they are comparable with
the dimensionless fluid rows; converged when max|residual| <= rtol.
"""
function BorderedFB(a, b0, G5, rows, u0; x0 = zeros(size(a, 1)),
                    rtol = 1e-10, maxit = 100, α = 0.5, β = 1e-4, τmin = 1e-12)
    N = size(a, 1); m = length(u0); Mt = N + m
    x = copy(x0); u = copy(u0)
    sc = max(norm(b0, Inf), norm(G5 * u0, Inf), 1.0)
    y = similar(x); yk = similar(x); xk = similar(x)
    F = zeros(Mt); Fk = zeros(Mt)
    px = zeros(N); qy = zeros(N); J = zeros(Mt, Mt)

    function eval!(F, y, x, u)
        mul!(y, a, x); y .+= b0; mul!(y, G5, u, 1.0, 1.0)
        @inbounds for i in 1:N
            F[i] = (sqrt(x[i]^2 + y[i]^2) - x[i] - y[i]) / sc
        end
        r, Ry, Ru = rows(y, u)
        F[N+1:Mt] .= r
        return Ry, Ru
    end

    Ry, Ru = eval!(F, y, x, u)
    f = 0.5 * dot(F, F)
    for it in 0:maxit
        maximum(abs, F) <= rtol && return x, u, y, it, true
        it == maxit && break
        # generalised Jacobian (at x = y = 0 take px = qy = -1, a valid element)
        @inbounds for i in 1:N
            si = sqrt(x[i]^2 + y[i]^2)
            if si < 1e-14 * sc
                px[i] = -1.0; qy[i] = -1.0
            else
                px[i] = x[i] / si - 1.0; qy[i] = y[i] / si - 1.0
            end
        end
        J[1:N, 1:N] .= (qy ./ sc) .* a
        @inbounds for i in 1:N
            J[i, i] += px[i] / sc
        end
        J[1:N, N+1:Mt] .= (qy ./ sc) .* G5
        J[N+1:Mt, 1:N] .= Ry * a
        J[N+1:Mt, N+1:Mt] .= Ry * G5 .+ Ru
        d = try
            J \ (-F)
        catch e
            e isa SingularException || rethrow()
            qr(J, ColumnNorm()) \ (-F)       # singular Jacobian (seen when the tail first closes)
        end
        τ = 1.0
        while true
            @views xk .= max.(x .+ τ .* d[1:N], 0.0)
            uk = u .+ τ .* d[N+1:Mt]
            Ryk, Ruk = eval!(Fk, yk, xk, uk)
            fk = 0.5 * dot(Fk, Fk)
            if fk <= (1 - 2β * τ) * f || τ < τmin
                x .= xk; u .= uk; y .= yk; F .= Fk; Ry, Ru = Ryk, Ruk; f = fk
                break
            end
            τ *= α
        end
    end
    return x, u, y, maxit, false
end

"""
    (Dn, Dss, Dds, Fluid) = FluidCrackSolve(P1,P2,P3,ν,G,λ,MidPoint,FaceNormalVector,HSFlag,
                                            MixedBC, FractureFlag, V0, ρ0, K; pa, pshift, g)

Shared core of the FluidVolume / FluidMass slip calculators. Friction on closed elements µ = 0.1,
Sf = 0 (as before). Fluid = (P = pressure at zref per crack, ρ̄ = mean fluid density per crack,
iterations, converged).
"""
function FluidCrackSolve(P1, P2, P3, ν, G, λ, MidPoint, FaceNormalVector, HSFlag, MixedBC, FractureFlag,
                         V0, ρ0, K; pa = 0.0, pshift = 0.0, g = 9.81)
    n = length(FaceNormalVector[:, 1])
    FixedEls = zeros(size(FractureFlag))
    # influence matrix and tractions from the remote stress / applied tractions
    (Dn, Dss, Dds, A, b) = SlipCalculator3D(P1, P2, P3, ν, G, λ, MidPoint, FaceNormalVector, HSFlag, MixedBC, FixedEls)
    (Area, HalfPerimeter) = AreaOfTriangle3D(P1[:, 1], P1[:, 2], P1[:, 3], P2[:, 1], P2[:, 2], P2[:, 3], P3[:, 1], P3[:, 2], P3[:, 3])
    Area = vec(Area)
    Scl = mean(HalfPerimeter) / G                 # scaling improves the solver
    A = A .* Scl
    Ainv = inv(A)
    # friction LCP (Kaven et al. 2012), µ and Sf as previously used for fluid-filled cracks
    µ = fill(0.1, n); Sf = zeros(n)
    (FricMat, FricVec, L1, L2, L3, L4, L5, ScX, D0) =
        SlipCalculator3D(Scl, n, InfMat(Ainv), BoundaryConditionsVec(b), µ, Sf)[1:9]
    a = FricMat.A
    b0 = copy(FricVec.b); b0[1:3n] .+= D0
    # fluid columns, equations, start
    cracks = MakeFluidCracks(FractureFlag, MidPoint, V0, ρ0, K)
    Gm = FluidColumns(Ainv, MidPoint, cracks, n; g = g)
    G5 = zeros(5n, size(Gm, 2)); G5[1:3n, :] .= Gm
    rows = FluidRows(cracks, Area, MidPoint, Scl, 5n; g = g, pshift = pshift, pa = pa)
    u0 = FluidLinearGuess(D0[1:n], Gm[1:n, :], cracks, Area, Scl)
    (x, u, y, its, ok) = BorderedFB(a, b0, G5, rows, u0)
    ok || @warn "Fluid-filled crack solve did not converge in $its iterations"
    Dn = y[L1] .* Scl
    Dss = (y[L2] .- x[L4]) .* Scl
    Dds = (y[L3] .- x[L5]) .* Scl
    nf = length(cracks)
    Fluid = (P = u[1:nf], ρ̄ = u[nf+1:2nf] .* [c.ρ0 for c in cracks], iterations = its, converged = ok)
    return Dn, Dss, Dds, Fluid
end
