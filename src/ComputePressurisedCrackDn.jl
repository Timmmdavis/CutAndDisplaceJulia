#If we have alreay found the optimum
function ComputePressurisedCrackDn(x::Tractions,Flag,B_old,Ainv,Scl,Area,Pcalc,n,Volume,
    ReturnVol,NumOfFractures,FricMatPrepped,FricVectorWithoutDisp,L1,L2,L3,L4,L5,
    D,Vects,Arrys,Flts,Ints,Mats,Bls,IntArrys)

#Setting up for the log descent (opening of parts of surface) 
#x - constant fluid pressure defined by the simulated annealing algo
#Logical - flag of cracks (0 is not a crack, 1:n are cracks). 
#B - traction boundary conditions at each element in the model
#Ainv - inverted traction influence matrix
#Scl - Scales the results (helps convergence of the friction solver). 
#Area - Area of each Triangles face. 
#PCalc - Pressure (max normal traction closing any element in model, scales
#the input x value accordingly).
#Volume - volumes of each crack (list). 

#@info x Flag B_old n Ainv[1:10,:]

#disp(x) # to diplay the current value

x=x.Tn

#Making sure no replacement happens inside func
B=copy(B_old);

#Add pressure to each seperate crack (different) 
for i=1:NumOfFractures #For each crack
    
    #Get the Current crack
    Indx=Flag.==i;
    Indx=findall(Indx)
    tmp=x[i]*Pcalc[i];
    B[Indx].+=tmp;  
    
end   

#Equation
D=Ainv*B;

##Extract arrays
Dn_=D[1:n];
Dss=D[n+1:2*n];
Dds=D[n*2+1:3*n];

# Way one (can have friction but ~2s per loop for 300+ tris)
if any(Flag.==0)
    println("it looks like you have elements that represent topography, deal with these correctly before passing into friction solver (see Davis 2019 ppr)")
end
if any(Dn_.>0)

    FricVectorWithoutDispIn=BoundaryConditionsVec(copy(FricVectorWithoutDisp.b))
    (Vects,Arrys,Flts,Ints,Mats,Bls,IntArrys)=FischerNewton.ResetInitArrays(Vects,Arrys,Flts,Ints,Mats,Bls,IntArrys);
    (Dn,Dss,Dds)=SlipCalculator3D(FricMatPrepped,FricVectorWithoutDispIn,L1,L2,L3,L4,L5,Scl,D,
                                  Vects,Arrys,Flts,Ints,Mats,Bls,IntArrys)

else
    println("Ding ding ding")
    Dn=Dn_.*0; #Already closed (wont open, we dont care)
    Dss=Dss.*0;
    Dds=Dds.*0;
end
#SCALES INSIDE!

#=
## Way two (same result as one, no friction involved)
# If we have negative volumes force these disps to negative
for i=eachindex(NumOfFractures) #For each crack
    #Get the volume of the current crack
    Vol_i_Expected=Volume[i];    
    if Vol_i_Expected<0
        println("DnFlipped")
        Dn_[Flag==i]=-Dn_[Flag==i];
    end
end

if any(Dn_.>0)
    #Get interpenetrating tris 
    Interpen=Dn_.<0; 
    #Parts of TnDn in Ainv that will be put to 0
    InterpenIndx=findall(Interpen);


    println("Trying to make it")
    I=1
    while any(Dn_.<0)
        I=I+1
        for j=1:length(InterpenIndx)

            indx=InterpenIndx[j]

            for i=1:n
                #Drop cols to 0 so closed elements opening have no opening influence on others opening (can still slip).
                #  A (  row  |  col ) = 0
                AinvF[indx,i]=0;
                AinvF[indx+n,i]=0;
                AinvF[indx+2*n,i]=0;
            end
            
            #Drop rows to 0 so closed elements cannot open (can still slip).
            #  A (  row  |  col ) = 0
            #AinvF[indx,:]=AinvF[indx,:].*0;

        end

        #Recompute D
        #D=Ainv*B;
        mul!(D,AinvF,B) 

        #Get results
        Dn_=D[1:n];
        Dss=D[n+1:2*n];
        Dds=D[n*2+1:3*n];
        if I > 20 #Max Iters
            break
        end
    end
    Dn=Dn_;
else
    Dn=Dn_.*0; #Already closed (wont open, we dont care)
end

#Scaling back the data. 
Dn=Dn.*Scl;
Dss=Dss.*Scl;
Dds=Dds.*Scl;
=#

return Dn,Dss,Dds

end


function ComputePressurisedCrackDn(x,Flag,B,Ainv,Scl,Area,Pcalc,n,Volume,
    ReturnVol,NumOfFractures,FricMatPrepped,FricVectorWithoutDisp,L1,L2,L3,L4,L5,
    D,Vects,Arrys,Flts,Ints,Mats,Bls,IntArrys)

CurrentPressure=Tractions(x,[],[]);
(Dn,Dss,Dds)=ComputePressurisedCrackDn(CurrentPressure,Flag,B,Ainv,Scl,Area,Pcalc,n,Volume,
    ReturnVol,NumOfFractures,FricMatPrepped,FricVectorWithoutDisp,L1,L2,L3,L4,L5,
    D,Vects,Arrys,Flts,Ints,Mats,Bls,IntArrys)

#Init some parameters before loop. 
X=0.; 
#Now compute objective function checking each crack matches given value
for i=1:NumOfFractures #For each crack
    
    #Get the Current crack
    Indx=findall(Flag.==i);
    #Get the volume of the current crack
    Vol_i_Expected=Volume[i];
    #Compute volume
    Vol_i_Computed=sum(Dn[Indx].*Area[Indx]);
    if ReturnVol==1
        X=Vol_i_Computed;
    else
        #Now compute the objective function 
        X=X+abs(Vol_i_Expected-Vol_i_Computed)^2;
    end
end

return X

end

