function R_3 = rot_3(theta)
% DESCRIPTION
%   Compute the rotation matrix for a rotation about the z-axis by an angle 
%   theta. This function is used in coordinate transformations involving a 
%   rotation around the third axis (z-axis).
%
% INPUTS        size    Type
%   theta       (1,1)   Double  Angle of rotation [radians]
%
% OUTPUTS       size    Type
%   R_3         (3,3)   Double  Rotation matrix for rotation about z-axis
%
% FUNCTION
%
    c = cos(theta);
    s = sin(theta);
    
    R_3 = [c  s  0;
          -s  c  0;
           0  0  1];
end