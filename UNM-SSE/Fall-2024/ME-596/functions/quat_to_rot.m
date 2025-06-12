function R = quat_to_rot(q)
    % DESCRIPTION
    %   Convert a quaternion to a rotation matrix.
    %
    % INPUTS       size    Type
    %   q          (4,1)   Double  Quaternion in the form [qx, qy, qz, qw]'
    %                        where qw is the scalar part, and qx, qy, qz 
    %                        are the vector components.
    %
    % OUTPUTS      size    Type
    %   R          (3,3)   Double  Rotation matrix
    %
    % FUNCTION
    
    q_4 = q(4);
    q = q(1:3);

    R = (q_4^2 - q' * q) * eye(3) + 2 * q * q' - 2 * q_4 * skew(q);
end
