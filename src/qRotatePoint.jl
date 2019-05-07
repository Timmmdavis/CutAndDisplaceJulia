function qRotatePoint( P, Qrotation )
# qRotatePoint: rotate a point according to rotation quaternion
# Protated = qRotatePoint( P, Qrotation )
# IN: 
#     P - point which is rotated
#     Qrotation - quaternion describing the rotation (axis and angle)
# 
# OUT:
#     Protated - rotated point 
#
# EXAMPLE:
#     Rotate point (1;2;3) around vector (4;5;6) by an angle of pi/2
#     P = [1;2;3];  # create the point
#     V = [4;5;6];  # create vector around which rotation is performed
#     Qrot = qGetRotQuaternion( pi/2, V );
#     P2 = qRotatePoint( P, Qrotate );  
#     
# VERSION: 03.03.2012

P = reshape( P, 3, 1 );
Q1 = [ 0; P ];
Q = qMul( Qrotation, Q1)
Q = qMul( Q, qInv(Qrotation) );
Protated = Q[2:4];

return Protated
end

function qMul( Q1, Q2 )
# qMul: quaternion multiplication 
# IN: 
#     Q1 - first quaternion
#     Q2 - second quaternion
# 
# OUT:
#     Q - output quaternion, Q = Q1*Q2
#     
# REMARKS:
#     1) Quaternion multiplication is not commutative, i.e. Q1*Q2 != Q2*Q1
#     2) Quaternion multiplication is associative, i.e. Q1*Q2*Q3 = Q1*(Q2*Q3)=(Q1*Q2)*Q3
# 
# VERSION: 03.03.2012

s1 = Q1[1];
s2 = Q2[1];
v1 = Q1[2:4];
v2 = Q2[2:4];

s =s1*s2 - dot( v1,v2);
v = s1*v2 + s2*v1 + cross( v1, v2 );
v = reshape( v, 3, 1 );
Q = [s;v];
return Q 
end   

function qInv( Q1 )
# qInv: quaternion reciprocal (inverse)
# Q = qInv( Q1 )
# IN: 
#     Q1 - input quaternion
# 
# OUT:
#     Q - reciprocal of Q1, i.e. Q1*Q = 1
#     
# VERSION: 03.03.2012

Q1 = reshape( Q1, 4, 1 );
Q = [Q1[1]; -Q1[2:4]];

Q = Q ./ ( sqrt( sum( Q1 .* Q1 )); )^2;

return Q
end