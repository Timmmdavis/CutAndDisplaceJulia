#Stress intensity around an elliptical crack under internal pressure P0 plus a linear pressure
#gradient (P1). Copied from VelocityEllipticalCrackAscent (StressIntensityAndOpeningEllipticalCrack.jl)
#so the elliptical-gradient test does not need that package.
#
#Formulas from Atroshchenko et al. 2009:
#Stress intensity factor for an embedded elliptical crack under arbitrary normal loading.
#Atroshchenko, University of Waterloo 2010 (thesis): Stress Intensity Factors for Elliptical and
#Semi-Elliptical cracks subjected to an arbitrary mode I loading.
#
#xc, zc: semi-axes. Assign_c1c2_andEllipticIntegrals -> shape terms; CompP0P1 -> the P0 and P1·zc
#that give KI = 0 at the base and Kc at the top; Comp_KC_EllipseEdge -> KI (pressure part,
#gradient part) at angle phi round the edge (from the z-axis).

function Assign_c1c2_andEllipticIntegrals(xc,zc)

        if zc>=xc
            c1=zc
            c2=xc
        elseif xc>zc
            c1=xc
            c2=zc
        end

        #Recalc for new crack shape
        k=sqrt((c1^2)-(c2^2))/c1;
        k_prime=c2/c1;
        (Kk,Ek) = Elliptic.ellipke(k^2);

        if zc>xc
            Q=1-2*k^2
        elseif zc<xc
            Q=1+k^2
        elseif xc==zc #Circular crack...
            Q=(pi + 2)/pi #Such that the term: (Q*Ek-Kk*(k_prime^2)) == 1
            k=sqrt(4/(3*pi)) #So it reduces to K around circular crack under grad
        end

    return c1,c2,k,k_prime,Kk,Ek,Q

end

function Comp_KC_EllipseEdge(c1,c2,zc,xc,Ek,Kk,k,k_prime,Q,P0,P1c,phi)
    #P0 internal pressure
    #P1c = gradient (P1) * zhalflength (zc)

    #angle away from z-axis
    alpha=phi

    if zc>xc
        P1c=-P1c
    else
        phi=phi+pi/2 #angle away from Atronshenko x-axis (fig 3.1/3.2) - we want from z-axis
    end

    #Where on the ellipses edge...
    angles=(((c2^2)/(c1^2))*(cos(phi)^2)+(sin(phi)^2))^(1/4);

    KI_InternalPressure=((sqrt(pi*c2)*P0)/Ek)*angles;

    KI_GradientPressure=((sqrt(pi*c2)*P1c*(k^2)*cos(alpha))/(Q*Ek-Kk*(k_prime^2)))*angles   

    return KI_InternalPressure,KI_GradientPressure

end

function CompP0P1(c1,c2,zc,xc,Ek,Kk,k,k_prime,Q,Kc)
    #Conditions: 0 at base, Kc at upper tip...
    #P1c = gradient (P1) * zhalflength (zc)

    if xc>=zc
        P0= (Kc*Ek)/(2*sqrt(pi*c2))
        P1c= (Kc*(Q*Ek-Kk*(k_prime^2))*c1)/(2*sqrt(pi*zc)*xc*(k^2))
    else
        P0= (Kc*Ek*sqrt(c1))/(2*sqrt(pi)*c2)
        P1c=-(Kc*(Q*Ek-Kk*(k_prime^2))*c1)/(2*sqrt(pi*zc)*xc*(k^2))
    end

    #Technically P1_zc that comes out... (P1 gradient * zc)
    return P0,P1c

end
