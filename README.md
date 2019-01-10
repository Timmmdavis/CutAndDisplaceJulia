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
| 2500  | 50*50   |186.39  | 2.534 |73.6 |
| 100   | 2*2     |  0.057 | 0.002 |28.5 |
| 100   | 10*10   |  0.086 | 0.004 |21.5 |
| 100   | 50*50   |  0.907 | 0.104 | 8.7 |
| 100   | 100*100 |  3.155 | 0.413 | 7.6 |
| 4     | 2*2     |   0.003056 |   0.00072 | 4.2  | 
| 16    | 4*4     |   0.013156 |   0.002   | 6.6  | 
| 100   | 10*10   |   0.088804 |   0.05    | 1.8  |
| 400   | 20*20   |   1.303736 |   0.72    | 1.8  |
| 900   | 30*30   |  10.439721 |   4.13    | 2.5  |
| 2500  | 50*50   | 197.946119 |  29.08    | 6.8  |
| 3600  | 60*60   | 723.903385 |  58.81    | 12.3 |
| 6400  | 80*80   |2906.401573 | 177.39    | 16.4 |

Relative speedup being the times from MATLAB(secs)/Julia(secs).
TD.jl produces all 3 inf matricies so is technically doing 3 times more (3x faster than reported above).  
