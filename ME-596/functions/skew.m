function A = skew(v)
% DESCRIPTION
%   Compute the skew-symmetric matrix of a 3-element vector. This function 
%   takes a vector and returns its corresponding skew-symmetric matrix, 
%   which is often used in cross-product calculations and rotational dynamics.
%
% INPUTS        size    Type
%   v           (3,1)   Double  Input 3-element vector
%
% OUTPUTS       size    Type
%   A           (3,3)   Double  Skew-symmetric matrix corresponding to the input vector
%
% FUNCTION
%
    % Check if the input is a 3-element vector
    if length(v) ~= 3
        error('Input must be a 3-element vector.');
    end
    
    % Create the skew-symmetric matrix
    A = [0, -v(3), v(2);
         v(3), 0, -v(1);
         -v(2), v(1), 0];
end