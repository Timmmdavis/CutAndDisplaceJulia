#If we have alreay found the optimum
function ComputePressurisedCrackDn(x::Tractions,Flag,B_old,Ainv,Scl,Area,Pcalc,n,Volume,ReturnVol,NumOfFractures)

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


#disp(x) # to diplay the current value

x=x.Tn
#Making sure no replacement happens inside func
B=copy(B_old);

#Add pressure to each seperate crack (different) 
for i=1:NumOfFractures #For each crack
    
    #Get the Current crack
    Indx=Flag.==i;

    for j=eachindex(Indx)
        if Indx[j]==true
            B[j]=B[j]+(x[i]*Pcalc[i]);  #Frak1 
       end
    end
    
end   

#Equation
D=Ainv*B;

#Extract arrays
Dn_=D[1:n];
Dss=D[n+1:2*n];
Dds=D[n*2+1:3*n];



#=# # Way one (can have friction but ~2s per loop for 300+ tris)
if any(Flag.==0)
    println("it looks like you have elements that represent topography, deal with these correctly before passing into friction solver (see Davis 2019 ppr)")
end
if any(Dn_.>0)
    Ainv=InfMat(Ainv); 
    b=BoundaryConditionsVec(B);
    µ=zeros(n);println("Setting arbitary fric params")
    Sf=zeros(n);
    (Dn,Dss,Dds)=SlipCalculator3D(Scl,n,Ainv,b,µ,Sf)
else
    Dn=Dn_.*0; #Already closed (wont open, we dont care)
    Dss=Dss.*0;
    Dds=Dds.*0;
end
#SCALES INSIDE!=#


## Way two (same result as one, no friction involved)
# If we have negative volumes force these disps to negative
for i=eachindex(NumOfFractures) #For each crack
    #Get the volume of the current crack
    Vol_i_Expected=Volume[i];    
    if Vol_i_Expected.<0
        println("DnFlipped")
        Dn_[Flag==i]=-Dn_[Flag==i];
    end
end

if any(Dn_.>0)
    #Get interpenetrating tris 
    Interpen=Dn_.<0; 
    #Parts of TnDn in Ainv that will be put to 0
    InterpenIndx=findall(Interpen);

    for i=1:n
        for j=eachindex(Interpen)
            if Interpen[j]==true
                #Drop cols to 0 so closed elements opening have no opening influence on others opening (can still slip).
                #  A (  row  |  col ) = 0
                Ainv[i,j]=0;
                Ainv[i,j+n]=0;
                Ainv[i,j+2*n]=0;
                #Drop rows to 0 so closed elements cannot open (can still slip).
                #  A (  row  |  col ) = 0
                Ainv[j,:]=Ainv[j,:].*0;
            end
        end
    end
    #Recompute D
    D=Ainv*B;
    #Get results
    Dn=D[1:n];
    Dss=D[n+1:2*n];
    Dds=D[n*2+1:3*n];
else
    Dn=Dn_.*0; #Already closed (wont open, we dont care)
end

#Scaling back the data. 
Dn=Dn*Scl;
Dss=Dss*Scl;
Dds=Dds*Scl;

return Dn,Dss,Dds

end


function ComputePressurisedCrackDn(x,Flag,B,Ainv,Scl,Area,Pcalc,NUM,Volume,ReturnVol,NumOfFractures)

CurrentPressure=Tractions(x,[],[]);
(Dn,Dss,Dds)=ComputePressurisedCrackDn(CurrentPressure,Flag,B,Ainv,Scl,Area,Pcalc,NUM,Volume,ReturnVol,NumOfFractures)

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
        X=X+abs(Vol_i_Expected-Vol_i_Computed);
    end
end

return X

end

