function [B] = Hessian_fun (x)
%  [B] = Hessian_fun (x)
%
% Example:
%   x = [1; 2; 3; 4; 5; 6];
%   B = Hessian_fun (x)
%   B - B.'

B = zeros(length(x), length(x));

% Odd indices:
for i=1:2:length(x)-1
  B(i,i)    = 2 - 40*(x(i+1) - x(i)^2) - 40*x(i)*(-2)*x(i);
  B(i+1,i)  = -40*x(i);
end    

% Even indices:
for i=2:2:length(x)
 B(i,i)    = 20;
 B(i-1, i) = -40*x(i-1);
end    

