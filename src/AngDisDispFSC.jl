function AngDisDispFSC(y1,y2,y3,cosB,sinB,cotB,cotB2,b1,b2,b3,nu,a)
# AngDisDispFSC calculates the harmonic function contribution to the 
# displacements associated with an angular dislocation in an elastic 
# half-space.

#sinB = sin(beta);
#cosB = cos(beta);
#cotB = cot(beta);
y3b = y3 +2*a;
z1b = y1*cosB+y3b*sinB;
z3b = -y1*sinB+y3b*cosB;
r2b = y1 ^2 +y2 ^2 +y3b^2;
rb = sqrt(r2b);

Fib = 2*atan(-y2 /(-(rb+y3b)*cotB2+y1)); # The Burgers' function

v1cb1 = b1/4/pi/(1 -nu)*(-2*(1 -nu)*(1 -2*nu)*Fib*cotB^2 +(1 -2*nu)*y2 /
    (rb+y3b)*((1 -2*nu-a/rb)*cotB-y1 /(rb+y3b)*(nu+a/rb))+(1 -2*nu)*
    y2 *cosB*cotB/(rb+z3b)*(cosB+a/rb)+a*y2 *(y3b-a)*cotB/rb^3 +y2 *
    (y3b-a)/(rb*(rb+y3b))*(-(1 -2*nu)*cotB+y1 /(rb+y3b)*(2*nu+a/rb)+
    a*y1 /rb^2)+y2 *(y3b-a)/(rb*(rb+z3b))*(cosB/(rb+z3b)*((rb*
    cosB+y3b)*((1 -2*nu)*cosB-a/rb)*cotB+2*(1 -nu)*(rb*sinB-y1)*cosB)-
    a*y3b*cosB*cotB/rb^2));

v2cb1 = b1/4/pi/(1 -nu)*((1 -2*nu)*((2*(1 -nu)*cotB^2 -nu)*log(rb+y3b)-(2*
    (1 -nu)*cotB^2 +1 -2*nu)*cosB*log(rb+z3b))-(1 -2*nu)/(rb+y3b)*(y1*
    cotB*(1 -2*nu-a/rb)+nu*y3b-a+y2 ^2 /(rb+y3b)*(nu+a/rb))-(1 -2*
    nu)*z1b*cotB/(rb+z3b)*(cosB+a/rb)-a*y1 *(y3b-a)*cotB/rb^3 +
    (y3b-a)/(rb+y3b)*(-2*nu+1 /rb*((1 -2*nu)*y1*cotB-a)+y2 ^2 /(rb*
    (rb+y3b))*(2*nu+a/rb)+a*y2 ^2 /rb^3)+(y3b-a)/(rb+z3b)*(cosB^2 -
    1 /rb*((1 -2*nu)*z1b*cotB+a*cosB)+a*y3b*z1b*cotB/rb^3 -1 /(rb*
    (rb+z3b))*(y2 ^2*cosB^2 -a*z1b*cotB/rb*(rb*cosB+y3b))));

v3cb1 = b1/4/pi/(1 -nu)*(2*(1 -nu)*(((1 -2*nu)*Fib*cotB)+(y2 /(rb+y3b)*(2*
    nu+a/rb))-(y2*cosB/(rb+z3b)*(cosB+a/rb)))+y2 *(y3b-a)/rb*(2*
    nu/(rb+y3b)+a/rb^2)+y2 *(y3b-a)*cosB/(rb*(rb+z3b))*(1 -2*nu-
    (rb*cosB+y3b)/(rb+z3b)*(cosB+a/rb)-a*y3b/rb^2));

v1cb2 = b2/4/pi/(1 -nu)*((1 -2*nu)*((2*(1 -nu)*cotB^2 +nu)*log(rb+y3b)-(2*
    (1 -nu)*cotB^2 +1)*cosB*log(rb+z3b))+(1 -2*nu)/(rb+y3b)*(-(1 -2*nu)*
    y1*cotB+nu*y3b-a+a*y1*cotB/rb+y1 ^2 /(rb+y3b)*(nu+a/rb))-(1 -2*
    nu)*cotB/(rb+z3b)*(z1b*cosB-a*(rb*sinB-y1)/(rb*cosB))-a*y1 *
    (y3b-a)*cotB/rb^3 +(y3b-a)/(rb+y3b)*(2*nu+1 /rb*((1 -2*nu)*y1*
    cotB+a)-y1 ^2 /(rb*(rb+y3b))*(2*nu+a/rb)-a*y1 ^2 /rb^3)+(y3b-a)*
    cotB/(rb+z3b)*(-cosB*sinB+a*y1 *y3b/(rb^3*cosB)+(rb*sinB-y1)/
    rb*(2*(1 -nu)*cosB-(rb*cosB+y3b)/(rb+z3b)*(1 +a/(rb*cosB)))));
                
v2cb2 = b2/4/pi/(1 -nu)*(2*(1 -nu)*(1 -2*nu)*Fib*cotB^2 +(1 -2*nu)*y2 /
    (rb+y3b)*(-(1 -2*nu-a/rb)*cotB+y1 /(rb+y3b)*(nu+a/rb))-(1 -2*nu)*
    y2*cotB/(rb+z3b)*(1 +a/(rb*cosB))-a*y2 *(y3b-a)*cotB/rb^3 +y2 *
    (y3b-a)/(rb*(rb+y3b))*((1 -2*nu)*cotB-2*nu*y1 /(rb+y3b)-a*y1 /rb*
    (1 /rb+1 /(rb+y3b)))+y2 *(y3b-a)*cotB/(rb*(rb+z3b))*(-2*(1 -nu)*
    cosB+(rb*cosB+y3b)/(rb+z3b)*(1 +a/(rb*cosB))+a*y3b/(rb^2*cosB)));
                
v3cb2 = b2/4/pi/(1 -nu)*(-2*(1 -nu)*(1 -2*nu)*cotB*(log(rb+y3b)-cosB*
    log(rb+z3b))-2*(1 -nu)*y1 /(rb+y3b)*(2*nu+a/rb)+2*(1 -nu)*z1b/(rb+
    z3b)*(cosB+a/rb)+(y3b-a)/rb*((1 -2*nu)*cotB-2*nu*y1 /(rb+y3b)-a*
    y1 /rb^2)-(y3b-a)/(rb+z3b)*(cosB*sinB+(rb*cosB+y3b)*cotB/rb*
    (2*(1 -nu)*cosB-(rb*cosB+y3b)/(rb+z3b))+a/rb*(sinB-y3b*z1b/
    rb^2 -z1b*(rb*cosB+y3b)/(rb*(rb+z3b)))));

v1cb3 = b3/4/pi/(1 -nu)*((1 -2*nu)*(y2 /(rb+y3b)*(1 +a/rb)-y2*cosB/(rb+
    z3b)*(cosB+a/rb))-y2 *(y3b-a)/rb*(a/rb^2 +1 /(rb+y3b))+y2 *
    (y3b-a)*cosB/(rb*(rb+z3b))*((rb*cosB+y3b)/(rb+z3b)*(cosB+a/
    rb)+a*y3b/rb^2));
                
v2cb3 = b3/4/pi/(1 -nu)*((1 -2*nu)*(-sinB*log(rb+z3b)-y1 /(rb+y3b)*(1 +a/
    rb)+z1b/(rb+z3b)*(cosB+a/rb))+y1 *(y3b-a)/rb*(a/rb^2 +1 /(rb+
    y3b))-(y3b-a)/(rb+z3b)*(sinB*(cosB-a/rb)+z1b/rb*(1 +a*y3b/
    rb^2)-1 /(rb*(rb+z3b))*(y2 ^2*cosB*sinB-a*z1b/rb*(rb*cosB+y3b))));
                
v3cb3 = b3/4/pi/(1 -nu)*(2*(1 -nu)*Fib+2*(1 -nu)*(y2*sinB/(rb+z3b)*(cosB+
    a/rb))+y2 *(y3b-a)*sinB/(rb*(rb+z3b))*(1 +(rb*cosB+y3b)/(rb+
    z3b)*(cosB+a/rb)+a*y3b/rb^2));

v1 = v1cb1[1]+v1cb2[1]+v1cb3[1];
v2 = v2cb1[1]+v2cb2[1]+v2cb3[1];
v3 = v3cb1[1]+v3cb2[1]+v3cb3[1];

return(v1,v2,v3)
end