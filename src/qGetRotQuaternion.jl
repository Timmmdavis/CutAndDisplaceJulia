function qGetRotQuaternion( teta, rotationVector )
# qCreateRotQuaternion: outputs a quaternion which is used to perform
# rotation
# Q = qGetRotQuaternion( teta, rotationVector )
# IN: 
#     teta - rotation angle
#     rotationVector - vector around which the rotation will be performed
# 
# OUT:
#     Q - rotation quaternion
#     
# VERSION: 03.03.2012

norm = sqrt(sum( rotationVector .* rotationVector ));
if( length( rotationVector ) != 3 )
    println( "rotationVector should have 3 coordinates!" );
    return;
end

rotationVector = reshape( rotationVector, 3, 1 );
if( norm > 0 )
    v = rotationVector / norm;
    Q = [ cos( teta/2 ); v*sin( teta/2 )];
else
    println("rotationVector cannot be 0 ");
end
    
return Q
end