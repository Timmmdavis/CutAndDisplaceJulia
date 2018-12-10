function CalcSideVec(PA,PB)

# Calculate TD side vector and the angle of the angular dislocation pair
SideVec = PB-PA;
eZ = [0;0;1];

G=-SideVec'*eZ/norm(SideVec);
beta = acos(G[1]);

return(SideVec,eZ,beta)

end