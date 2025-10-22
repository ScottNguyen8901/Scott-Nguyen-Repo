function [alpha] = backtrack_line_search (fun, f_x_k, x_k, p_k, g, initial_alpha, c, rho)
%
% [alpha] = backtrack_line_search (fun, f_x_k, x_k, p_k, g, initial_alpha, c, rho)
% Computes a value for alpha so that:
%   f(x_k + alpha*p_k) < f(x_k) + c*alpha*grad(f)^T*p_k,
% This is accomplished by reducing alpha to alpha*rho in each iteration. 
%
% Input:
%   fun:           the function to be minimized.
%   f_x_k:         is the function evaluated at x_k.
%   initial_alpha: the initial value for alpha.
%                  This value for alpha should be as large as possible,
%                  so that a sufficient reduction is expected there.
%   c, rho:        minimization parameters (see above).
%   g:             the gradient of the function evaluated at x_k.
%
% Output:
%   alpha:        the step length value that satisfies the condition (see above).
%
% Examples:
%   fun    = @(x) x^2-1;
%   g_fun  = @(x) 2*x;
%   initial_alpha = 1;
%   c   = 0.1;
%   rho = 0.9
%   x_k = -1;
%   p_k = 10;
%   f_x_k = fun(x_k);
%   g_k   = g_fun(x_k);
%   [alpha] = backtrack_line_search (fun, f_x_k, x_k, p_k, g_k, initial_alpha, c, rho); 
%  x = -3:0.01:3;
%  plot(x, x.^2-1), hold on;
%  plot(x_k, fun(x_k), '+');
%  plot(x_k + alpha*p_k, fun(x_k + alpha*p_k), '*');
%  

Max_iterations = 100;

% Procedure 3.1 from Nocedal and Wright, page 41.
k = 1;
alpha = initial_alpha;
while ((fun(x_k + alpha*p_k) > f_x_k+c*alpha*(g.')*p_k) && ... % sufficient decrease?
       (k<Max_iterations))
  alpha = rho*alpha;
  k = k+1;
end

