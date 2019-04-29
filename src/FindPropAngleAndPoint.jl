function FindPropAngleAndPoint( Fe,FaceNormalVector,I,MaxStrainEnergy,AvgTriangleEdgeLength )
#UNTITLED4 Summary of this function goes here
#   Detailed explanation goes here

    FeMd=Fe.FeMd[I,:];
    FeEv=Fe.FeEv[I,:];
    FeM2Ev=Fe.FeM2Ev[I,:];
    #FeLe=Fe.FeLe[I,:];
    StrainEnergy=Fe.StrainEnergy[I,:]
    K1=Fe.K1[I,:];K1=K1[1];
    K2=Fe.K2[I,:];K2=K2[1];
    NormalV=FaceNormalVector[I,:]

    L=AvgTriangleEdgeLength#FeLe*sind(60); #Side length of tri - Max propagation
    ScaledL=(L/MaxStrainEnergy)*StrainEnergy;

    #Could be cleverer here, "C:\Users\timmm\Documents\PapersPDF\igabem3d_01doubleSpace_StressIntensity.pdf"
    #Eq.43 (paris law)

    #Right lateral movement so look in neg theta domain (see PF P375)
    if K2[1]>0.0 
        θo = range(-0.5*pi,stop=0.0,length=90);
        #Eq 9.78 Pollard:
        σθθ=(K1*sin.(θo)).+(K2*((3.0*cos.(θo)).-1.0));
        v=σθθ; x=θo;
        #using Interpolations
        itp = interpolate((v,),x, Gridded(Linear()))
        #finding σθθ at θo=0
        Ang=TryInterp(itp,0.0)


    #Left lateral movement
    else
        θo = range(0.0,stop=0.5*pi,length=90);
        #Eq 9.78 Pollard:
        σθθ=(K1*sin.(θo)).+(K2*((3.0*cos.(θo)).-1.0));
        v=σθθ; x=θo;
        #using Interpolations
        itp = interpolate((v,),x, Gridded(Linear()))
        #finding σθθ at θo=0
        Ang=TryInterp(itp,0.0)   


    end
    #From igabem3d_01doubleSpace_StressIntensity
    #θc=2.0*atan((-2.0*(K2/K1))/(1.0+sqrt(1.0+8.0*(K2/K1)^2.0))) 
    #Doesnt appear to work...
    #@info Ang θc


    #New Pnt - crack local coords
    X=1.0; #along crack
    Y=0.0; #crack normal dir
    Z=0.0; #crack front edge dir
    

    #@info Ang
    (X,Y)=CutAndDisplaceJulia.RotateObject2D!(X,Y,0.0,0.0,cos(Ang),sin(Ang))
    #Now we have a vector length 1 pointing in right direction
    #local crack coords

    (X,Y,Z) = RotateObject3DNewCoords!(X,Y,Z,0.0,0.0,0.0,FeM2Ev,NormalV,FeEv);
    #FeM2EvK=normr([X Y Z]);
    FeM2EvK=[X Y Z]

    #Rotated and length 1 into world coords

    #Check if two vectors match direction:
    Vect=(dot(FeM2Ev',FeM2EvK'))<=0;
    #Flip if not
    if Vect==1
        FeM2EvK.=-FeM2EvK; 
    end
    
    
    #Scale by length away from tip in cartesian coords
    DistFromMd=FeM2EvK.*ScaledL;
    #So new point location is:
    NwPntCX=FeMd[1].+DistFromMd[1] 
    NwPntCY=FeMd[2].+DistFromMd[2] 
    NwPntCZ=FeMd[3].+DistFromMd[3];


return NwPntCX,NwPntCY,NwPntCZ,Ang
end

