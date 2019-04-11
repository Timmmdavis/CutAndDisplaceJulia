function FindPropAngleAndPoint( FeMd,FeM2Ev,FeLe,FeEv,NormalV,K2,K1 )
#UNTITLED4 Summary of this function goes here
#   Detailed explanation goes here

    L=FeLe*sind(60); #Side length of tri
    #Could be cleverer here, "C:\Users\timmm\Documents\PapersPDF\igabem3d_01doubleSpace_StressIntensity.pdf"
    #Eq.43 (paris law)

    #Right lateral movement so look in neg theta domain (see PF P375)
    if K2>0 
        ThetaZ = range(-0.5*pi,stop=0.0,length=90);
        #Eq 9.78 Pollard:
        Null=(K1.*sin.(ThetaZ))+(K2.*((3.0.*cos.(ThetaZ)).-1.0));
        x=Null; v=ThetaZ;
        #using Interpolations
        itp = interpolate((x,),v, Gridded(Linear()))
        ThetaZ=TryInterp(itp,0.0)

    #Left lateral movement
    else
        ThetaZ = range(0.0,stop=0.5*pi,length=90);
        #Eq 9.78 Pollard:
        Null=(K1.*sin.(ThetaZ))+(K2.*((3.0.*cos.(ThetaZ)).-1.0));
        x=Null; v=ThetaZ;
        #using Interpolations
        itp = interpolate((x,),v, Gridded(Linear()))
        ThetaZ=TryInterp(itp,0.0)                   
    end

    #Flip EV if normal is flipped:
    FeM2EvN=FeM2Ev;
    Bad=findall(NormalV[3]<0);

    FeEv[Bad,:]=-FeEv[Bad,:];
    FeM2EvN[Bad,:]=-FeM2EvN[Bad,:];

    #New Pnt= 
    X=1.0; 
    Z=0.0;
    Y=tan(ThetaZ);
    (X,Y,Z) = RotateObject3DNewCoords!(X,Y,Z,0.0,0.0,0.0,FeM2Ev,NormalV,FeEv);
    FeM2EvK=normr([X,Y,Z]);
           
    #Check if two vectors match direction:
    Vect=(dot(FeM2Ev',FeM2EvK'))<=0;
    #Flip if not
    if Vect==1
        FeM2EvK=-FeM2EvK; 
    end

    #Length in cartesian coords
    DistFromMd=FeM2EvK.*L;
    #So new point location
    NwPntCX=FeMd[1]+DistFromMd[1] 
    NwPntCY=FeMd[2]+DistFromMd[2] 
    NwPntCZ=FeMd[3]+DistFromMd[3];


return NwPntCX,NwPntCY,NwPntCZ,ThetaZ
end
