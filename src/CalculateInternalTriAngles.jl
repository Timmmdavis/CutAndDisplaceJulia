function CalculateInternalTriAngles(P1,P2,P3)
#returns internal triangle angles in degrees

n=size(P1,1);

#Prepping for loop
IntAngA=zeros(n);
IntAngB=zeros(n);
IntAngC=zeros(n);

P1toP2Vec=normr(P2.-P1)
P2toP1Vec=normr(P1.-P2)

P1toP3Vec=normr(P3.-P1)
P3toP1Vec=normr(P1.-P3)

P2toP3Vec=normr(P3.-P2)
P3toP2Vec=normr(P2.-P3)

for i = 1:n

    IntAngA[i]=acos(dot(P1toP2Vec[i,:],P1toP3Vec[i,:]))
    IntAngB[i]=acos(dot(P2toP1Vec[i,:],P2toP3Vec[i,:]))
    IntAngC[i]=acos(dot(P3toP1Vec[i,:],P3toP2Vec[i,:]))

end

IntAngA=rad2deg.(IntAngA)
IntAngB=rad2deg.(IntAngB)
IntAngC=rad2deg.(IntAngC)

return IntAngA,IntAngB,IntAngC

end


