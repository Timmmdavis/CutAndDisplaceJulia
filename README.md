# CutAndDisplaceInJulia

[![Build Status](https://travis-ci.com/Timmmdavis/CutAndDisplaceJulia.svg?token=1HhESyMNyqzV8R22Pqq6&branch=master)](https://travis-ci.com/Timmmdavis/CutAndDisplaceJulia)
[![Coverage Status](https://codecov.io/gh/Timmmdavis/CutAndDisplaceJulia/branch/master/graph/badge.svg?token=IbbZ8n4385)](https://codecov.io/gh/Timmmdavis/CutAndDisplaceJulia)

This repository contains CutAndDisplace working in Julia. 
The test folder contains some contains some simple benchmarks against analytical solutions in hope of catching regressions (run through Travis on each commit). Running through the full benchmarks in https://github.com/Timmmdavis/CutAndDisplace would be desirable to catch issues related to coordinates, sign etc. 

In src the main functions are 
LD.jl 
and 
TD.jl
These compute stresses and displacements due to analytical line and triangular dislocations. 

LD.jl I havent tested for unstable types etc, TD.jl I have. If hooked up to code it would be good to do thos. 
For more improvements in TD.jl the functions CalculateLocalTriCoords and GlobalToTDECoords could be improved in terms of memory management. Also Parts in the 'FS' functions that are enabled when the 'ImageFlag' is set to 1. 
I recommend removing Threads.@threads from the outer loop of these functions before debugging/working on the code.  

## Halfspace speed comparison (Triangular dislocations)

For commit 25dffe8 in terms of speed for TD.jl compared to MATLAB (CompareHalfSpaceSol.m) using funcs from: https://github.com/Timmmdavis/CutAndDisplace/blob/master/3dCode/TDFunctions/TDstrainHS.m commit b07067a  (HS strain):
(Windows machine - Intel(R) Xeon(R) CPU X5472 @ 3.00GHz	Base speed:	2.99 GHz 	Sockets:	2	Cores:	8	Logical processors:	8). Using Btime (BenchmarkTools): 

| No of Tris    | No of obs Points | MATLAB (seconds)  | Julia (seconds) | Relative speedup |
| ------------- |:----------------:| -----------------:| -------------:  | --------------:  |
| 2     | 50*50   |  0.37  |  0.06  |  5.9 |
| 10    | 50*50   |  1.69  |  0.12  | 13.4 |
| 100   | 50*50   | 17.12  |  0.87  | 19.6 |
| 2500  | 50*50   |425.94  | 21.38  | 19.9 |
| 100   | 2*2     |  0.72  |  0.003 | 240     |
| 100   | 10*10   |  4.09  |  0.03  | 120     |
| 100   | 50*50   | 17.12  |  0.87  | 19.6 |
| 100   | 100*100 | 59.84  |  3.51  | 17 |

## Fullspace speed comparison (Triangular dislocations)

In terms of speed for TD.jl against FS MATLAB (CompareHalfSpaceSol.m with HS=0;)

| No of Tris    | No of obs Points | MATLAB (seconds)  | Julia (seconds) | Relative speedup |
| ------------- |:----------------:| -----------------:| -------------:  | --------------:  |
| 2     | 50*50   |  0.013 | 0.005 | 2.6 | 
| 10    | 50*50   |  0.072 | 0.013 | 5.5 | 
| 100   | 50*50   |  0.907 | 0.104 | 8.7 |
| 2500  | 50*50   |186.39  | 4.12 | 45 |
| 100   | 2*2     |  0.057 | 0.005 |11.4|
| 100   | 10*10   |  0.086 | 0.04 | 2.15 |
| 100   | 50*50   |  0.907 | 0.104 | 8.7 |
| 100   | 100*100 |  3.155 | 0.413 | 7.6 |
| 4     | 2*2     |   0.003 |  0.0001307 | 23  | 
| 16    | 4*4     |   0.013 |  0.00049   | 27  | 
| 100   | 10*10   |   0.09 |   0.00722   | 12 |
| 400   | 20*20   |   1.30 |   0.1031    | 13  |
| 900   | 30*30   |  10.44 |   0.518     | 20  |
| 2500  | 50*50   | 197.95 |  4.12       | 48  |
| 3600  | 60*60   | 723.90 |  8.64      | 84 |
| 6400  | 80*80   |2906.40 | 27.58      | 105 |


Relative speedup being the times from MATLAB(secs)/Julia(secs).
TD.jl produces all 3 inf matricies so is technically doing 3 times more (3x faster than reported above).  


## Installing PyCall etc windows:

Install Julia 1.0.3 64 bit (https://julialang.org/downloads/)
Make sure you can call from cmd.exe, i.e. 
```
>julia.exe 
```
runs it. Do this by setting the julia.exe dir as a path enviroment variable. 

Install Python Windows 64-x86 (https://www.python.org/downloads/release/python-372/). 
Again make sure its callable from cmd, i.e. 
```
> python 
```

Open julia and call 
```
] add PyCall
] build PyCall
```

Now in cmd call - 
```
> python -m pip install julia
```

In cmd call python
```
> python
```
In Python call 
```
>>> import julia
>>> j = julia.Julia()
```
If you get an error here it promts you to where a new python.exe is installed. 
For me this was -> C:\Users\timmm\.julia\conda\3
Change your python enviroment variable to match this and check you can call python from cmd, do prompts it recommends. 
Then repeat steps for python from above again. 

Now:
In cmd call python
```
> python
```
In Python call 
```
>>> import julia
>>> j = julia.Julia()
>>> from julia import Base
>>> Base.sind(90)
>>> from julia import TravisTest #https://github.com/Timmmdavis/JuliaTravisTest
>>> TravisTest.MrsFunc(2)
```
1.0


## MATLAB 2 Python 2 Julia
Making sure that the python.exe that works with julia (conda) on the cmd line is the one that is your env var

Script.py =
```
import sys
import julia
import numpy as np
import json
import scipy.io

#Show we can print in MATLAB command window
#print("helloworld")

#See what we have brought in
#print(dir())

#This is the bit where we catch the input vars... 
if __name__ == '__main__':
    x = np.array(json.loads(sys.argv[1]))
    y = np.array(json.loads(sys.argv[2]))
    #z = float(sys.argv[3]) #If you have floats etc...

#Check the inputs have been picked up correctly
#print(locals())

#Compute something in Julia (from a module)
j = julia.Julia()
j.include("WorkOnArrays.jl")
z=j.WorkOnArrays(x,y)

#Write results in MATLAB command window
sys.stdout.write(str(z))

#Or save as .Mat file
scipy.io.savemat('PythonOutput.mat', dict(x=x, y=y, z=z))
```

then in MATLAB using the system cmd call (In the dir of the script)   
```
>> clear
>> a=[[0,0,0,0];[0,1,1,0];[0,1,1,1];[0,0,0,0]]
>> b=[[0,0,1,0];[0,3,1,0];[0,1,1,1];[0,0,0,0]]
>> a=Mat2NumpyStr(a);
>> b=Mat2NumpyStr(b);
>> [status,cmdout]=system(['python Script.py '... 
 a ...
 ' '...
 b])
and output will be your result (as a string)
Alternativly load the matrix
>> load('PythonOutput.mat');
where Mat2NumpyStr.m:
```
function [str]=Mat2NumpyStr(mat)  

str=num2str(mat,'%.1f,');
for i=1:size(str,1)
    str(i,end)=';';
end
str = reshape(str.',1,[]);
str=['[[' str ']]'];
str = strrep(str,';','],[');
str = strrep(str,',[]','');
str = strrep(str,' ','');
```
and Script.py
```
import sys
import julia
import numpy as np
import json
import scipy.io

#Show we can print in MATLAB command window
#print("helloworld")

#See what we have brought in
#print(dir())

#This is the bit where we catch the input vars... 
if __name__ == '__main__':
    x = np.array(json.loads(sys.argv[1]))
    y = np.array(json.loads(sys.argv[2]))
    #z = float(sys.argv[3]) #If you have floats etc...

#Check the inputs have been picked up correctly
#print(locals())

#Compute something in Julia (from a module)
j = julia.Julia()
j.include("WorkOnArrays.jl") #in the current folder in MATLAB
z=j.WorkOnArrays(x,y)

#Write results in MATLAB command window
sys.stdout.write(str(z))

#Or save as .Mat file
scipy.io.savemat('PythonOutput.mat', dict(x=x, y=y, z=z))
```
and 
WorkOnArrays.jl
```
function WorkOnArrays(X::Array,Y::Array)
Z=zeros(size(X));
for i=1:length(X)
	#Rotate to new axes Ax Ay Az
	Z[i]=(X[i]-rand())+(Y[i]-rand());
end
return(Z)
end
```

## TD example script 
```
import sys
import julia
import numpy as np
import json
import scipy.io

#See what we have brought in
#print(dir())

#This is the bit where we catch the input vars... 
if __name__ == '__main__':
	x = np.array(json.loads(sys.argv[1]))
	y = np.array(json.loads(sys.argv[2]))
	z = np.array(json.loads(sys.argv[3]))
	P1 = np.array(json.loads(sys.argv[4]))
	P2 = np.array(json.loads(sys.argv[5]))
	P3 = np.array(json.loads(sys.argv[6]))
	DssVec = np.array(json.loads(sys.argv[7]))
	DdsVec = np.array(json.loads(sys.argv[8]))
	DnVec =  np.array(json.loads(sys.argv[9]))	
	nu=float(sys.argv[10])
	mu=float(sys.argv[11])
	DispFlag=int(sys.argv[12])
	StressFlag=int(sys.argv[13])
	HSFlag=int(sys.argv[14])
	
#Check the inputs have been picked up correctly
#print(locals())

#Compute something in Julia (from a module)
j = julia.Julia()
from julia import MyModule
(ExxDn,EyyDn,EzzDn,ExyDn,ExzDn,EyzDn,
ExxDss,EyyDss,EzzDss,ExyDss,ExzDss,EyzDss,
ExxDds,EyyDds,EzzDds,ExyDds,ExzDds,EyzDds,
UxDn,UyDn,UzDn,UxDss,UyDss,UzDss,UxDds,UyDds,UzDds)=MyModule.TD(x,y,z,P1,P2,P3,DssVec,DdsVec,DnVec,nu,mu,DispFlag,StressFlag,HSFlag)

#Or save as .Mat file
scipy.io.savemat('PythonOutput.mat', dict(
ExxDn=ExxDn,EyyDn=EyyDn,EzzDn=EzzDn,ExyDn=ExyDn,ExzDn=ExzDn,EyzDn=EyzDn,
ExxDss=ExxDss,EyyDss=EyyDss,EzzDss=EzzDss,ExyDss=ExyDss,ExzDss=ExzDss,EyzDss=EyzDss,
ExxDds=ExxDds,EyyDds=EyyDds,EzzDds=EzzDds,ExyDds=ExyDds,ExzDds=ExzDds,EyzDds=EyzDds,
UxDn=UxDn,UyDn=UyDn,UzDn=UzDn,UxDss=UxDss,UyDss=UyDss,UzDss=UzDss,UxDds=UxDds,UyDds=UyDds,UzDds=UzDds))
```

## TD Setup script 

```
x = linspace(-10,10,2);
[x,y]=meshgrid(x,x);
z=ones(size(x))*-0; %Ground surface
x=reshape(x,numel(x),1);
y=reshape(y,numel(y),1);
z=reshape(z,numel(z),1);

X=[-1  1 -1 -1 1 1];
Y=[ 1 -1 -1 -1 1 1];
Z=[ 0  -3  0  -3 0 0];
P1=[X(1) Y(1) Z(1);X(2) Y(2) Z(2)];
P2=[X(5) Y(5) Z(5);X(4) Y(4) Z(4)];
P3=[X(3) Y(3) Z(3);X(6) Y(6) Z(6)];

DdsVec=[2.0,2.0];
DnVec=[0.6,0.6];
DssVec=[0.5,0.5];

mu=1.; %const 
nu=0.25; %const 

DispFlag=1;
StressFlag=1;
HSFlag=0;


x=Mat2NumpyStr(x);
y=Mat2NumpyStr(y);
z=Mat2NumpyStr(z);

P1=Mat2NumpyStr(P1);
P2=Mat2NumpyStr(P2);
P3=Mat2NumpyStr(P3);

DdsVec=Mat2NumpyStr(DdsVec);
DnVec=Mat2NumpyStr(DnVec);
DssVec=Mat2NumpyStr(DssVec);

nu=num2str(nu,'%.1f');
mu=num2str(mu,'%.1f');

DispFlag=num2str(DispFlag);
StressFlag=num2str(StressFlag);
HSFlag=num2str(HSFlag);


[status,cmdout]=system(['python Script_TD.py '... 
 x ...
 ' '...
 y ...
 ' '... 
 z ...
 ' '...
 P1 ...
 ' '...  
 P2 ...
 ' '...  
 P3 ...
 ' '...  
 DdsVec ...
 ' '...  
 DnVec ...
 ' '...  
 DssVec ...
 ' '...  
 nu ...
 ' '...  
 mu ...
 ' '...   
 DispFlag ...
 ' '...  
 StressFlag ...
 ' '...    
 HSFlag])
```