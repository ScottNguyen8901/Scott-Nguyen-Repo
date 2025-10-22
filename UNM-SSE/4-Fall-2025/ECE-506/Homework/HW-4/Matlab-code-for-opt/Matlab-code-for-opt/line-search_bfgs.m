function [x, x_all] = BFGS (fun, G_fun, initial_x, initial_H, vec_eps)
%
% [x, x_all] = BFGS (fun, G_fun, initial_x, initial_H, vec_eps)
% Minimizes the given function fun using the BFGS algorithm,
% a quasi-Newton method.
%
% Input:
%    fun:        function to optimize.  
%    G_fun:      a function for estimating the gradient at a point.
%    initial_x:  initial value for x.
%    initial_H:  initial value for H.
%    vec_eps:    gradient vector magnitude value for terminating.
%
% Output:
%    x:     the computed value where the function is minimized.
%    x_all: the trajectory of all the points that were taken.
%
% Example:
%  fun   = @(x) x(1).^2 + 5*x(2).^2;
%  g_fun = @(x) [2*x(1); 10*x(2)];
%  initial_x = [-2; -3];
%  vec_eps   = 0.01;
%  initial_H = [2 0; 0 10];
%  [x, x_all] = BFGS (fun, g_fun, initial_x, initial_H, vec_eps);
%  [x1, x2] = meshgrid(-4:0.01:4, -4:0.01:4);
%  [cs,h]=contour(x1, x2, x1.^2 + 5*x2.^2); hold on;
%  plot(x_all(1,:), x_all(2,:));
%  plot(x_all(1,:), x_all(2,:), '*');
  

% page 198

% Initialize line-search parameters:
initial_a_i = 1;
a_max = 20;
c1    = 1.0e-4;
c2    = 0.9;

% Initialize the setup:
x_k = initial_x;
H_k = initial_H;
g_k = G_fun(x_k);

% Initialize the iterations. 
Max_iterations = 100;
k = 1;

N     = length(x_k);
Ident = diag(ones(1,N));

% Continue until the gradient has been sufficiently reduced:
x_all = x_k;
while (vec_norm(g_k)>vec_eps) && (k<Max_iterations)
 
  p_k = - H_k * g_k;
  
  % Build phi, phi':
  f_k     = fun(x_k);
  phi     = @(a) fun(x_k + a*p_k);
  phi_der = @(a) ((G_fun(x_k + a*p_k)).')*p_k; % See (A.16), page 583.
  [a_opt] = line_search (phi, phi_der, f_k, ((g_k).')*p_k, initial_a_i, a_max, c1, c2);
  x_k_plus_1 = x_k + a_opt*p_k;
  x_all = [x_all x_k_plus_1];
  
  % Compute the update parameters:
  g_k_plus_1 = G_fun(x_k_plus_1);
  s_k = x_k_plus_1 - x_k;
  y_k = g_k_plus_1 - g_k;  

  % Verify that y_k^T*s_k>0
  rho_k = 1 / ((y_k.')*s_k);
  if (rho_k <= 0)
      disp('y_k^T*s_k < 0 => exiting ...');
      break
  end
  
  % Update the inverse:
  A   = Ident - rho_k*s_k*(y_k.');
  A_t = Ident - rho_k*y_k*(s_k.');
  H_k_plus_1 = A*H_k*A_t + rho_k*s_k*(s_k.');

  % Update the old parameters here:
  x_k = x_k_plus_1;
  g_k = g_k_plus_1;
  H_k = H_k_plus_1;
  k   = k +1;
end
x = x_k;