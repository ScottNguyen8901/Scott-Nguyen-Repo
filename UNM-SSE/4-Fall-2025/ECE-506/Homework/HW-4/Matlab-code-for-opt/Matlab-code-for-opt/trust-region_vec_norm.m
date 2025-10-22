function [val] = vec_norm (p)
%    [val] = vec_norm (p)
% Input:
%   p:   a vector of any size
% Output:
%   val: the norm of the vector.
% Example:
%   p = [1; 2];
%   disp(sprintf('Error = %f', vec_norm(p) - sqrt(5)));
%

val = sqrt((p.')*p);