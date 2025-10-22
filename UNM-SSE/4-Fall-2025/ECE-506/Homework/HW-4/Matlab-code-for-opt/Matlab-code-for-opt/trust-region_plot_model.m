function [] = plot_model (fun, G_fun, H_fun, model_x, fun_vec, x_vals, y_vals)
%
% [] = plot_model (fun, G_fun, H_fun, model_x, fun_vec, x_vals, y_vals)
% provides contour plots of fun and a quadratic model at model_x.
%
% Input:
%  fun:     a 2-D function for which we are building a model.
%  G_fun:   a gradient function that returns the gradient.
%  H_fun:   a Hessian function for fun.
%  fun_vec: a vectorized version of fun.
%  model_x: the value of x where we want to build the model.
%  x_vals, y_vals: contain [min_val, max_val] that specify where
%                  to plot the model.
% Output:
%  A contour plot of fun and the local quadratic model at model_x.
%
% Example:
%  % Provide scalar and vector version of the function:
%  fun     = @(x)     10*(x(2) - x(1)^2)^2 + (1-x(1))^2;
%  fun_vec = @(x1,x2) 10*(x2   - x1.^2).^2 + (1-x1).^2;
%  % Evaluate the gradient in G_fun and the Hessian in H_fun:
%  G_fun = @(x) [20*(x(2)-x(1)^2)*(-2*x(1)) + 2*(1-x(1))*(-1); ...
%                20*(x(2) - x(1)^2)];
%  H_fun = @(x) [-40*x(2)+120*x(1)^2+2, -40*x(1); -40*x(1) 20];
%  plot_model (fun, G_fun, H_fun, [0; -1], fun_vec, [-5 5], [-5 10]);
%

% Build the local quadratic model:
f_k = fun(model_x);
g_k = G_fun(model_x);
B_k = H_fun(model_x);

% Create a function centered at model_x:
model_fun = @(p) f_k + (g_k.')*p + 0.5*(p.')*B_k*p;

% Evaluate the model for a range of values:
N = 100;
min_x = x_vals(1);
max_x = x_vals(2);
min_y = y_vals(1);
max_y = y_vals(2);

Delta_x = (max_x - min_x)/(N-1);
Delta_y = (max_y - min_y)/(N-1);

x = zeros(N, N);
y = x;
z = x;
for i=1:N
 for j=1:N
     x(i,j) = min_x+(j-1)*Delta_x;
     y(i,j) = min_y+(i-1)*Delta_y;
     p = [x(i,j); y(i,j)];
     model_m(i,j) = model_fun(p);
 end
end
[x1,x2] = meshgrid(min_x:0.1:max_x, min_y:0.1:max_y);
contour(x1, x2, fun_vec(x1,x2), 100), hold on;
contour(x, y, model_m, 100)
plot(model_x(1), model_x(2), '*');