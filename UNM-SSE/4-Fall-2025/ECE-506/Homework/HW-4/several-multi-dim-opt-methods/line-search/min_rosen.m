close all
clear all

% Define Rosenbrock function:
%  1. Scalar version: R_fun, and
%  2. Vector version: R_vfun.
R_fun  = @(x)      100*(x(2) - x(1).^2).^2 + (1-x(1)).^2;
R_vfun = @(x1, x2) 100*(x2   - x1.^2).^2   + (1-x1).^2;

% Define the gradient and Hessian functions:
R_g_fun = @(x) [-400*x(1)*(x(2)-x(1)^2)-2+2*x(1); ...
                 200*(x(2)-x(1)^2)];
R_H_fun = @(x) [ -400*(x(2)-x(1)^2)+800*x(1)^2 + 2, ...
                 -400*x(1); -400*x(1), 200];
             
% Define minimum gradient for stopping the iterations:
vec_eps   = 0.01;

% Run the methods at an easy point:
initial_x = [1.2; 1.2];
initial_H = R_H_fun (initial_x);
[x_steep,  x_all_steep]  = steepest_descent (R_fun, R_g_fun, initial_x, vec_eps);
[x_newton, x_all_newton] = newton (R_fun, R_g_fun, R_H_fun, initial_x, vec_eps);
[x_BFGS,   x_all_BFGS]   = BFGS (R_fun, R_g_fun, initial_x, initial_H, vec_eps);


% Plot all the trajectories
figure(1);
x_range = 0.5:0.01:1.5;
y_range = x_range;
subplot(1,3,1);
plot_trajectory ('Steepest Descent', R_vfun, x_range, y_range, x_all_steep); 
subplot(1,3,2);
plot_trajectory ('Newton', R_vfun, x_range, y_range, x_all_newton);
subplot(1,3,3);
plot_trajectory ('BFGS',   R_vfun, x_range, y_range, x_all_BFGS);


% Run the methods at an easy point:
initial_x = [-1.2; 1];
initial_H = R_H_fun (initial_x);
[x_steep,  x_all_steep]  = steepest_descent (R_fun, R_g_fun, initial_x, vec_eps);
[x_newton, x_all_newton] = newton (R_fun, R_g_fun, R_H_fun, initial_x, vec_eps);
[x_BFGS,   x_all_BFGS]   = BFGS (R_fun, R_g_fun, initial_x, initial_H, vec_eps);


% Plot all the trajectories
figure(2);
x_range = -1.5:0.01:1.5;
y_range = -0.5:0.01:1.5;
subplot(1,3,1);
plot_trajectory ('Steepest Descent', R_vfun, x_range, y_range, x_all_steep); 
subplot(1,3,2);
plot_trajectory ('Newton', R_vfun, x_range, y_range, x_all_newton);
subplot(1,3,3);
plot_trajectory ('BFGS',   R_vfun, x_range, y_range, x_all_BFGS);
