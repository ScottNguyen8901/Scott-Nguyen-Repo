clear all
close all

% Provide scalar and vector version of the function:
fun     = @(x)     10*(x(2) - x(1)^2)^2 + (1-x(1))^2;
fun_vec = @(x1,x2) 10*(x2   - x1.^2).^2 + (1-x1).^2;

% Evaluate the gradient in G_fun and the Hessian in H_fun:
G_fun = @(x) [20*(x(2)-x(1)^2)*(-2*x(1)) + 2*(1-x(1))*(-1); ...
              20*(x(2) - x(1)^2)];
H_fun = @(x) [-40*x(2)+120*x(1)^2+2, -40*x(1); -40*x(1) 20];

% Plot the model at [0; -1]:
figure(1);
initial_x = [0; -1];
plot_model (fun, G_fun, H_fun, initial_x, fun_vec, [-5 5], [-5 10]);

% Minimize using the dogleg method:
g = G_fun(initial_x);
B = H_fun(initial_x);
for Delta = 0.1:0.01:2
  [p_dogleg] = dogleg   (g, B, Delta);
  [min_p]    = min_quad (g, B, Delta);
  
  plot([initial_x(1) initial_x(1)+p_dogleg(1)], ...
       [initial_x(2) initial_x(2)+p_dogleg(2)]);
  plot([initial_x(1) initial_x(1)+min_p(1)], ...
       [initial_x(2) initial_x(2)+min_p(2)], ':');
end    

% Repeat at [0; 0.5]:
figure(2);
initial_x = [0; 0.5];
plot_model (fun, G_fun, H_fun, initial_x, fun_vec, [-5 5], [-5 10]);

% Minimize using the dogleg method:
g = G_fun(initial_x);
B = H_fun(initial_x);
for Delta = 0.1:0.01:2
  [p_dogleg] = dogleg (g, B, Delta);
  [min_p] = min_quad (g, B, Delta);
  
  plot([initial_x(1) initial_x(1)+p_dogleg(1)], ...
       [initial_x(2) initial_x(2)+p_dogleg(2)]);
  plot([initial_x(1) initial_x(1)+min_p(1)], ...
       [initial_x(2) initial_x(2)+min_p(2)], ':');
end    
