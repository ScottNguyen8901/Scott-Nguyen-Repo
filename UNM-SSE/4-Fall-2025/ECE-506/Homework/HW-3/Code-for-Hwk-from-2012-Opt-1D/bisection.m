function [x_vals, fun_evals] = bisection(a, b, tols, fun, max_iters)
%
% [x_vals, fun_evals] = bisection(a, b, tols, fun, max_iters)
%
% Input:
%   a, b:      specify an interval for the Bisection method.
%   tols:      tolerances: |a-b|>tols(1) and f(a),f(b)>tols(2).
%   fun:       the function to evaluate.
%   max_iters: maximum allowed iterations. 
% Output:
%   x_vals:    an array holding all the generate values for x.
%   fun_evals: the number of function evaluations used in each iterations.
%
% Descrition:  
%   The program applies the bisection method and terminates if:
%     (sign(f(a)*f(b))<0) or (|a-b|<tol1) or (f(a),f(b)<tol2) or (iterations=max_iters)
% Example:
%   fun_sqr1 = @(x) x^2 - 3;
%   [a, b, max_iters] = deal(1.0, 2.1, 100);
%   tols = [1.0e-5 1.0e-8];
%   [x_vals, fun_evals] = bisection(a, b, tols, fun_sqr1, max_iters);
%   plot(fun_evals, abs(x_vals-sqrt(3)));
%

% initialize variables:
fa = feval(fun, a);
fb = feval(fun, b);
iterations  = 0;
x_vals      = [a];
fun_evals   = [2];
total_evals = 2;

tol1 = tols(1);
tol2 = tols(2);

while ((iterations<max_iters)& ...
       (abs(a-b)>tol1) & ...
       (abs(fa)>tol2)  & ...
       (abs(fb)>tol2)  & ...
       (sign(fa*fb)<0))
    
  mid_point = (a+b)/2;
  fun_mid   = feval(fun, mid_point);
 
  if (sign(fun_mid) ~= sign(fa))
     b  = mid_point;
     fb = fun_mid;
  else % Move a to the mid-point:
     a  = mid_point;
     fa = fun_mid;
  end
  
  total_evals = total_evals+1;
  x_vals      = [x_vals mid_point]; 
  fun_evals   = [fun_evals total_evals];
  iterations  = iterations+1;    
end    