function  TensorTransformation2D(Pxx,Pyy,Pxy,CosAx,CosAy )
#See 2D MATLAB func for more info. Is it faster with loop like for 3D?

#First make a variable to speed stuff up:
CosAxPxx=CosAx.*Pxy;

#Just doing matrix multiplication (shown below) but with indexing. 
P11q=(CosAx.*Pxx)+(-CosAy.*Pxy);
P12q=CosAxPxx+(-CosAy.*Pyy);
P21q=(CosAy.*Pxx)+CosAxPxx;
P22q=(CosAy.*Pxy)+(CosAx.*Pyy);

P11=(CosAx.*P11q)+(-CosAy.*P12q);
P12=((CosAx.*P12q)+(-CosAy.*P22q)+(CosAy.*P11q)+(CosAx.*P21q))./2;
P22=(CosAy.*P21q)+(CosAx.*P22q);


return(P11,P22,P12)
end