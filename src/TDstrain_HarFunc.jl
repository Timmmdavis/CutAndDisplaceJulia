function TDstrain_HarFunc(X,Y,Z,P1,P2,P3,by,bz,bx,mu,lambda,nu)
# TDstrain_HarFunc calculates the harmonic function contribution to the
# strains and stresses associated with a triangular dislocation in a 
# half-space. The function cancels the surface normal tractions induced by 
# the main and image dislocations.


# Calculate unit strike, dip and normal to TD vectors: For a horizontal TD 
# as an exception, if the normal vector points upward, the strike and dip 
# vectors point Northward and Westward, whereas if the normal vector points
# downward, the strike and dip vectors point Southward and Westward, 
# respectively.
(Vnorm,Vstrike,Vdip)=CalculateLocalTriCoords(P1,P2,P3)

## Transform slip vector components from TDCS into EFCS
(bX,bY,bZ) = RotateObject3DNewCoords(bx,by,bz,0,0,0,Vnorm,Vstrike,Vdip);

# Calculate contribution of angular dislocation pair on each TD side 
(Exx1,Eyy1,Ezz1,Exy1,Exz1,Eyz1) = AngSetupStrainFSC(X,Y,Z,bX,bY,bZ,P1,P2,mu,lambda,nu); # P1P2
(Exx2,Eyy2,Ezz2,Exy2,Exz2,Eyz2) = AngSetupStrainFSC(X,Y,Z,bX,bY,bZ,P2,P3,mu,lambda,nu); # P2P3
(Exx3,Eyy3,Ezz3,Exy3,Exz3,Eyz3) = AngSetupStrainFSC(X,Y,Z,bX,bY,bZ,P3,P1,mu,lambda,nu); # P3P1

# Calculate total harmonic function contribution to strains and stresses
Exx=Exx1.+Exx2.+Exx3;
Eyy=Eyy1.+Eyy2.+Eyy3;
Ezz=Ezz1.+Ezz2.+Ezz3;
Exy=Exy1.+Exy2.+Exy3;
Exz=Exz1.+Exz2.+Exz3;
Eyz=Eyz1.+Eyz2.+Eyz3;

return(Exx,Eyy,Ezz,Exy,Exz,Eyz);
end