function [x_vals, fun_evals] = newton(x0, typx, tols, fun, fun_der, max_iters)
%
% [x_vals, fun_evals] = newton(x0, tol, fun, fun_der, max_iters)
%
% Input:
%   x0:        specify a starting point for the search.
%   typx:      typical value expected for the root (see Dennis and
%              Schnabel, page 27).
%   tols:      tolerances: |f(x)|>tols(1), |f'(x)|>tols(2),                       
%                          |x_n+1 - x_n|/max(typx, |x_n+1|))>tols(3).
%   fun:       the function to evaluate.
%   fun_der:   the function derivative.
%   max_iters: maximum allowed iterations. 
% Output:
%   x_vals:    an array holding all the generate values for x.
%   fun_evals: the number of function evaluations used in each iterations.
%
% Descrition:  
%   The program applies the Newton's method with exact derivatives until
%   the tolerances are violated or iterations=max_iters.
% Example:
%   fun     = @(x) x^4   - 12*x^3 + 47*x^2 - 60*x;
%   fun_der = @(x) 4*x^3 - 36*x^2 + 94*x   - 60;
%   tols = [1.0e-7 1.0e-7 1.0e-7];
%   [x0, max_iters] = deal(-3, 100);
%   [x_vals, fun_evals]  = newton(x0, typx, tols, fun, fun_der, max_iters);
%   plot(fun_evals, abs(x_vals-0), '+');
%

% initialize variables:
fx    = feval(fun, x0);
f_der = feval(fun_der, x0);
iterations  = 0;
x_vals      = [x0];
fun_evals   = [2];
total_evals = 2;

x_prev = x0-0.001; % for checking in while loop.
x = x0;

% Setup all the tolerances:
[tol1, tol2, tol3] = deal(tols(1), tols(2), tols(3));

while ((iterations<max_iters)& ...
       (abs(fx)>tol1) & ...
       (abs(f_der)>tol2)  & ...
       (abs(x-x_prev)/max(typx,abs(x)) > tol3))

  x_prev = x;
  x = x - fx/f_der; % Main equation in Newton's algorithm!   
  
  fx    = feval(fun, x);
  f_der = feval(fun_der, x);
  
  % Update statistics:       
  total_evals = total_evals+2;
  x_vals      = [x_vals x]; 
  fun_evals   = [fun_evals total_evals];
  iterations  = iterations+1;    
end    