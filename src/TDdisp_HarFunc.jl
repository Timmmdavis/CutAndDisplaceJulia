function TDdisp_HarFunc(X,Y,Z,P1,P2,P3,by,bz,bx,nu)
# TDdisp_HarFunc calculates the harmonic function contribution to the
# displacements associated with a triangular dislocation in a half-space.
# The function cancels the surface normal tractions induced by the main and
# image dislocations.

#bx = Ts; # Tensile-slip
#by = Ss; # Strike-slip
#bz = Ds; # Dip-slip

# Calculate unit strike, dip and normal to TD vectors: For a horizontal TD 
# as an exception, if the normal vector points upward, the strike and dip 
# vectors point Northward and Westward, whereas if the normal vector points
# downward, the strike and dip vectors point Southward and Westward, 
# respectively.
(Vnorm,Vstrike,Vdip)=CalculateLocalTriCoords(P1,P2,P3)

## Transform slip vector components from TDCS into EFCS
(bX,bY,bZ) = RotateObject3DNewCoords(bx,by,bz,0,0,0,Vnorm,Vstrike,Vdip);

# Calculate contribution of angular dislocation pair on each TD side 
(u1,v1,w1) = AngSetupDispFSC(X,Y,Z,bX,bY,bZ,P1,P2,nu); # Side P1P2
(u2,v2,w2) = AngSetupDispFSC(X,Y,Z,bX,bY,bZ,P2,P3,nu); # Side P2P3
(u3,v3,w3) = AngSetupDispFSC(X,Y,Z,bX,bY,bZ,P3,P1,nu); # Side P3P1

# Calculate total harmonic function contribution to displacements
ue = u1+u2+u3;
un = v1+v2+v3;
uv = w1+w2+w3;

return(ue,un,uv);
end