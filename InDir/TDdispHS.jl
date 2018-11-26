function TDdispHS(X,Y,Z,P1,P2,P3,by,bz,bx,nu)
# TDdispHS 
# Calculates displacements associated with a triangular dislocation in an 
# elastic half-space.
#
# TD: Triangular Dislocation
# EFCS: Earth-Fixed Coordinate System
# TDCS: Triangular Dislocation Coordinate System
# ADCS: Angular Dislocation Coordinate System
# 
# INPUTS
# X, Y and Z: 
# Coordinates of calculation points in EFCS (East, North, Up). X, Y and Z 
# must have the same size. 
#
# P1,P2 and P3:
# Coordinates of TD vertices in EFCS. 
# 
# Ss, Ds and Ts:
# TD slip vector components (Strike-slip, Dip-slip, Tensile-slip).
#
# nu:
# Poisson's ratio.
#
# OUTPUTS
# ue, un and uv:
# Calculated displacement vector components in EFCS. ue, un and uv have
# the same unit as Ss, Ds and Ts in the inputs.
# 
# 
# Example: Calculate and plot the first component of displacement vector 
# on a regular grid.
# 
# [X,Y,Z] = meshgrid(-3:.02:3,-3:.02:3,-5);
# [ue,un,uv] = TDdispHS(X,Y,Z,[-1 0 0],[1 -1 -1],[0 1.5 -2],-1,2,3,.25);
# h = surf(X,Y,reshape(ue,size(X)),'edgecolor','none');
# view(2)
# axis equal
# axis tight
# set(gcf,'renderer','painters')

# Reference journal article: 
# Nikkhoo M. and Walter T.R., 2015. Triangular dislocation: An analytical, 
# artefact-free solution. 
# Submitted to Geophysical Journal International 

# Copyright (c) 2014 Mehdi Nikkhoo
# 
# Permission is hereby granted, free of charge, to any person obtaining a 
# copy of this software and associated documentation files 
# (the "Software"), to deal in the Software without restriction, including 
# without limitation the rights to use, copy, modify, merge, publish, 
# distribute, sublicense, and/or sell copies of the Software, and to permit
# persons to whom the Software is furnished to do so, subject to the 
# following conditions:
# 
# The above copyright notice and this permission notice shall be included 
# in all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
# NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
# USE OR OTHER DEALINGS IN THE SOFTWARE.

# I appreciate any comments or bug reports.

# Mehdi Nikkhoo
# created: 2013.1.24
# Last modified: 2014.7.30
# 
# VolcanoTectonics Research Group
# Section 2.1, Physics of Earthquakes and Volcanoes
# Department 2, Physics of the Earth
# Helmholtz Centre Potsdam
# German Research Centre for Geosciences (GFZ)
# 
# email: 
# mehdi.nikkhoo@gfz-potsdam.de 
# mehdi.nikkhoo@gmail.com


if any(Z .>0) | any(P1[3] .>0) | any(P2[3] .>0) | any(P3[3] .>0)
    error("Half-space solution: Z coordinates must be negative!")
end

#bx = Ts; # Tensile-slip
#by = Ss; # Strike-slip
#bz = Ds; # Dip-slip


#Init some vars
p1 = zeros(3,1);
p2 = zeros(3,1);
p3 = zeros(3,1);
P1_1=P1[1];P1_2=P1[2];P1_3=P1[3];
P2_1=P2[1];P2_2=P2[2];P2_3=P2[3];
P3_1=P3[1];P3_2=P3[2];P3_3=P3[3];

# Calculate main dislocation contribution to displacements
(ueMS,unMS,uvMS) = TDdispFS(X,Y,Z,P1,P2,P3,by,bz,bx,nu);

# Calculate harmonic function contribution to displacements
(ueFSC,unFSC,uvFSC) = TDdisp_HarFunc(X,Y,Z,P1,P2,P3,by,bz,bx,nu);

# Calculate image dislocation contribution to displacements
P1[3] = -P1[3];
P2[3] = -P2[3];
P3[3] = -P3[3];
(ueIS,unIS,uvIS) = TDdispFS(X,Y,Z,P1,P2,P3,by,bz,bx,nu);
if P1[3]==0 && P2[3]==0 && P3[3]==0
    uvIS = -uvIS;
end

# Calculate the complete displacement vector components in EFCS
ue = ueMS+ueIS+ueFSC;
un = unMS+unIS+unFSC;
uv = uvMS+uvIS+uvFSC;

if P1[3]==0 && P2[3]==0 && P3[3]==0
    ue = -ue;
    un = -un;
    uv = -uv;
end
DisplacementXYZ=[ue,un,uv];

return DisplacementXYZ
end