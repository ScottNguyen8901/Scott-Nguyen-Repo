function R = euler_to_rot(phi, a)
    % DESCRIPTION
    %   Compute the rotation matrix for a rotation by an angle phi around
    %   an arbitrary unit vector a. This function uses the Rodrigues' 
    %   rotation formula to compute the rotation matrix.
    %
    % INPUTS        size    Type
    %   phi        (1,1)   Double  Rotation angle in radians
    %   a          (3,1)   Double  Unit vector representing the axis of rotation
    %
    % OUTPUTS       size    Type
    %   R          (3,3)   Double  Rotation matrix
    %
    % FUNCTION

    R = eye(3) + sin(phi) * skew(a) + (1 - cos(phi)) * skew(a)^2;
end
