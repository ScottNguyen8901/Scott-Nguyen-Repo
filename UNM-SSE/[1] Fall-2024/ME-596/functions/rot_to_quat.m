function q = rot_to_quat(R)
    % DESCRIPTION
    %   Converts a rotation matrix to a quaternion representation.
    %
    % INPUTS       size    Type    Description
    %   R          (3,3)   Double  Rotation matrix
    %
    % OUTPUTS      size    Type    Description
    %   q          (4,1)   Double  Quaternion in the form [qx, qy, qz, qw]',
    %                              where qw is the scalar part and qx, qy, qz
    %                              are the vector components.
    %
    % FUNCTION

    % Compute the scalar part of the quaternion
    q_4 = (1 / 2) * sqrt(1 + trace(R));
    
    % Compute the vector part of the quaternion
    q = (1 / (4 * q_4)) * [(R(2,3) - R(3,2)); 
                           (R(3,1) - R(1,3)); 
                           (R(1,2) - R(2,1))];
    
    % Append the scalar part to the quaternion
    q = [q; q_4];
end
