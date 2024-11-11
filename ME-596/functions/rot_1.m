function R_1 = rot_1(theta)
    % DESCRIPTION
    %   Compute the rotation matrix for a rotation about the x-axis by an angle 
    %   theta. This function is used in coordinate transformations involving a 
    %   rotation around the first axis (x-axis).
    %
    % INPUTS        size    Type
    %   theta       (1,1)   Double  Angle of rotation [radians]
    %
    % OUTPUTS       size    Type
    %   R_1         (3,3)   Double  Rotation matrix for rotation about x-axis
    %
    % FUNCTION
    %
    c = cos(theta);
    s = sin(theta);
    
    R_1 = [1  0  0;
           0  c  s;
           0 -s  c];
end