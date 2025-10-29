function [x, x_all, counts] = newton (fun, g_fun, Hessian_fun, initial_x, vec_eps)
%
% [x, x_all] = newton (fun, g_fun, Hessian_fun, initial_x, vec_eps)
% Applies Newton's algorithm with line-search (proc. 3.1 of Nocedal & Wright)
% to minimize the function given in fun.
%
% Inputs:
%   initial_x:   an initial gues at the minimizing point.
%   fun:         the function to be minimized.
%   g_fun:       the gradient function, returning the gradient.
%   Hessian_fun: the Hessian function returning the Hessian.
%   vec_eps:     the threshold for the gradient magnitude. Below this value,
%                the gradient is considered to be zero.
% Outputs:
%   x:         returns the estimate of x that minimizes the function fun.
%   x_all:     returns all the estimates of x for minimizing the function
%              (for debugging).  
% Example:
% fun   = @(x) x(1).^2 + 5*x(2).^2;
% g_fun = @(x) [2*x(1); 10*x(2)];
% Hessian_fun = @(x) [2, 0; 0 10];
% initial_x = [-2; -3];
% vec_eps   = 0.01;
% [x, x_all] = newton (fun, g_fun, Hessian_fun, initial_x, vec_eps);
% [x1, x2] = meshgrid(-4:0.01:4, -4:0.01:4);
% [cs,h]=contour(x1, x2, x1.^2 + 5*x2.^2); hold on;
% %clabel(cs,h);
% plot(x_all(1,:), x_all(2,:));
% plot(x_all(1,:), x_all(2,:), '*');


% Local, steepest descent parameters:
initial_alpha = 1;
c   = 1e-4;
rho = 0.5;

max_iterations = 100;
k   = 1;

% Initialize and store the x
x_k   = initial_x;
x_all = x_k(:);

% --- added counters ---
counts = struct('f',0,'g',0,'H',0);

while (1)
   % Set p_k to the negative of the gradient of the function:
   g     = g_fun(x_k); counts.g = counts.g + 1;
   H     = Hessian_fun(x_k); counts.H = counts.H + 1;
   
   % NAIVE computation of p_k, only works for positive definite matrices!
   p_k   = -H\g;
   f_x_k = fun(x_k); counts.f = counts.f + 1;
   
   % Apply backtrach line search for minimizing the function.
   [alpha_k] = backtrack_line_search (fun, f_x_k, x_k, p_k, g, initial_alpha, c, rho);
   x_k = x_k + alpha_k*p_k;
   k = k+1;
  
  % Store the x:
  x_all = [x_all x_k(:)];
  
  % Print the step:
  % disp(sprintf('alpha_k = %f', alpha_k));

  % Stop the search for small gradient or a large number of iterations.
  if ((k == max_iterations) | (vec_norm(g)<vec_eps))
      x = x_k;
      break;
  end    
end