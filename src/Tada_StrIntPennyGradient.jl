function Tada_StrIntPennyGradient(Δρ,c,Theta)

#c is radius
#p in tada = Δρ*c when there is a gradient in stress that is Δρ

# Penny crack gradient - Tada Stress analysis handbook P.355:
K1=(4/(3*pi))*(Δρ*c)*sqrt(pi*c).*cos.(Theta)

return K1

end