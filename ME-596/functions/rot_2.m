function R_2 = rot_2(theta)
    %
    % DESCRIPTION
    %   Compute the rotation matrix for a rotation about the y-axis by an angle 
    %   theta. This function is commonly used in coordinate transformations 
    %   where a rotation around the second axis (y-axis) is needed.
    %
    % INPUTS        size    Type
    %   theta       (1,1)   Double  Angle of rotation [radians]
    %
    % OUTPUTS       size    Type
    %   R_2         (3,3)   Double  Rotation matrix for rotation about y-axis
    %
    % FUNCTION
    %

    c = cos(theta);
    s = sin(theta);
    
    R_2 = [c 0 -s;
           0 1 0;
           s 0 c];
end
