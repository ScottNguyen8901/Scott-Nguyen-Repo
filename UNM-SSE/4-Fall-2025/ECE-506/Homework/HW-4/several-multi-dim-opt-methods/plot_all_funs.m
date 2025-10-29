% Initialize environment.
clear all;
close all;

% Build a 2D coordinate system.
N1 = 10;
N2 = 10;
[x1, x2] = meshgrid(0:N1, 0:N2);

% Build three basic functions cases:
f_pos_def = @(x)   2*(x(1)-N1/2).^2 + (x(2)-N2/2).^2;
f_neg_def = @(x) - 2*(x(1)-N1/2).^2 - (x(2)-N2/2).^2;
f_mixed   = @(x)   2*(x(1)-N1/2).^2 - (x(2)-N2/2).^2;

% Visualize three basic functions:
num_of_lines = 20;
plot_fun_opt ('f pos. def.', f_pos_def, x1, x2, num_of_lines);
plot_fun_opt ('f neg.def.',  f_neg_def, x1, x2, num_of_lines);
plot_fun_opt ('f mixed',     f_mixed,   x1, x2, num_of_lines);


% Oscillatory behavior.
f_osc = @(x) 10*cos((2*2*pi/N1)*x(1))*cos((2*2*pi/N2)*x(2));
plot_fun_opt ('f osc.', f_osc, x1, x2, num_of_lines);


% Combine
f_comb = @(x) f_pos_def(x) + f_osc(x);
plot_fun_opt('f quad+osc', f_comb, x1, x2, num_of_lines);