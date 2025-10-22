function [x_vals, fun_evals] = quasi_newton2(x0, typx, tols, fun, max_iters)
%
% [x_vals, fun_evals] = quasi_newton2(x0, typx, tols, fun, max_iters)
%
% Input:
%   x0:        specify a starting point for the search.
%   typx:      typical value expected for the root (see Dennis and
%              Schnabel, page 27).
%   tols:      tolerances: |f(x)|>tols(1), |f'(x)|>tols(2),                       
%                          |x_n+1 - x_n|/max(typx, |x_n+1|))>tols(3).
%   fun:       the function to evaluate.
%   max_iters: maximum allowed iterations. 
% Output:
%   x_vals:    an array holding all the generate values for x.
%   fun_evals: the number of function evaluations used in each iterations.
%
% Descrition:  
%   The program applies the Quasi-Newton's method with finite-differencing,
%   where the derivative is approximated using:
%       f'(x) appox= (f(x_n) - f(x_n-1))/(x_n - x_n-1); % from the previous value.       
%   approximation to the derivatives until the tolerances are violated or
%   iterations=max_iters.
% Example:
%   fun = @(x) x^4   - 12*x^3 + 47*x^2 - 60*x;
%   [x0, max_iters] = deal(-3, 100);
%   tols = [1.0e-7 1.0e-7 1.0e-7];
%   typx = 1;
%   [x_vals, fun_evals]  = quasi_newton2(x0, typx, tols, fun, max_iters)
%   plot(fun_evals, abs(x_vals-0),'*');
%

% initialize variables:
x_prev = x0 - 1.0e-3;    % Previous value of x.
x      = x0;             % New x.

fx_prev = feval(fun, x_prev);
fx_new  = feval(fun, x);
f_der = (fx_new - fx_prev)/(x-x_prev);

iterations  = 0;
x_vals      = [x0];
fun_evals   = [2];
total_evals = 2;

% Setup all the tolerances:
[tol1, tol2, tol3] = deal(tols(1), tols(2), tols(3));

while ((iterations<max_iters)& ...
       (abs(fx_new)>tol1) & ...
       (abs(f_der)>tol2)  & ...
       (abs(x-x_prev)/max(typx,abs(x)) > tol3))
  % Store the previous iteration:  
  x_prev  = x;    
  fx_prev = fx_new;
  
  % Take the Newton step:
  x = x - fx_new/f_der; % Main equation in Newton's algorithm!   
  
  % Evaluate the new point:
  fx_new = feval(fun, x);
  f_der  = (fx_new - fx_prev)/(x - x_prev);  
  
  % Update statistics:       
  total_evals = total_evals+1;
  x_vals      = [x_vals x]; 
  fun_evals   = [fun_evals total_evals];
  iterations  = iterations+1;    
end    