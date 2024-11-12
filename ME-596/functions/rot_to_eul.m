function [phi, a] = rot_to_eul(R)
    % DESCRIPTION
    %   Compute the rotation angle and Euler axis from a given rotation matrix R.
    %   This function calculates the rotation angle phi (always positive) and the 
    %   unit vector a representing the Euler axis based on the properties of the rotation matrix.
    %
    % INPUTS        size    Type
    %   R           (3,3)   Double  Rotation matrix
    %
    % OUTPUTS       size    Type
    %   phi         (1,1)   Double  Rotation angle in radians (non-negative)
    %   a           (3,1)   Double  Unit vector representing the Euler axis (consistent direction)
    %
    % FUNCTION
    %   The function calculates:
    %     phi = acos((1/2) * (trace(R) - 1))
    %     a = 1 / (2 * sin(phi)) * (R' - R)
    
    % Calculate the rotation angle
    phi = acos((1/2) * (trace(R) - 1));
    
    % Calculate the skew-symmetric matrix for axis vector a
    a_skew = (1 / (2 * sin(phi))) * (R' - R);
    
    % Extract the axis vector a from the skew-symmetric matrix
    a = [a_skew(3,2); a_skew(1,3); a_skew(2,1)];

end
