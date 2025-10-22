function [] = test_zero_funs (fun, fun_der, the_roots, interval, N)
%
% test_zero_funs (fun, the_roots, interval, N)
%

% Generate random guesses for Bisection:
x_low = interval(1);
x_hi  = interval(2);
root1 = the_roots(1);

bis_a = root1 - rand([N 1])*(root1 - x_low) - 1.0e-6;
bis_b = root1 + rand([N 1])*(x_hi  - root1) + 1.0e-6;
max_iters = 1000;

% General Parameters for Newton-related algorithms:
tols = [1.0e-7 1.0e-7 1.0e-7]; % Tolerances for f(.), f'(.), and relative error.
typx = 1;                      % Typical value expected for optimal x.
delta_x = 1.0e-5;              % delta_x used for finite differencing.

for num = 1:N
  a = bis_a(num);
  b = bis_b(num);
  [bis_x, bis_evals] = bisection(a, b, tols, fun, max_iters);    

  % Use a for x0:
  x0 = a;
 [newt_x, newt_evals] = newton(x0, typx, tols, fun, fun_der, max_iters);  
 [q_newt_x,  q_newt_evals]  = quasi_newton(x0, typx, tols, fun, delta_x, max_iters);
 [q_newt2_x, q_newt2_evals] = quasi_newton2(x0, typx, tols, fun, max_iters);    
 [ds_x,      ds_evals]      = ds_method(x0, typx, tols, fun, max_iters);

  % Plot all the results
  plot_all_iters;
end    
