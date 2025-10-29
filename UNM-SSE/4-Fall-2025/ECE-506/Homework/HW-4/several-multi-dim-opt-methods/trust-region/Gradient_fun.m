function [g] = Gradient_fun (x)

% Expand-out the function:
%  (1 - x1)^2 + 10(x2 - x1^2)^2 +
%  (1 - x3)^2 + 10(x4 - x3^2)^2 +
%  (1 - x5)^2 + 10(x6 - x5^2)^2 +
%  (1 - x7)^2 + 10(x8 - x7^2)^2 +
%  ...
%  (1 - x_(2n-1))^2 + 10(x_(2n) - x_(2n-1)^2)^2
% Gives:
%  Odd-numbered ders:
%   -2(1-x_i) + 20(x_i+1 - x_i^2)*(-2x_i)
%  Even-numbered ders:
%   20(x_i - x_i-1^2)

% Odd-indexed components appear twice:  
for i=1:2:length(x)-1
 g(i) = -2*(1-x(i)) - 40*(x(i+1) - x(i)^2)*x(i); 
end    

% Even-indexed components appear once:
for i=2:2:length(x)
 g(i) = 20*(x(i) - x(i-1)^2);
end    
g = g(:);