% Regular multi-objective optimization methods:
%  Regular methods assume that objective values are available over
%  regularly spaced intervals. This is a trivial case and it is one for
%  which all of the problems can be solved directly.
clear
clc
run_reg_mo_obj

function run_reg_mo_obj
N = 10;
M = 10;
fs = 14; % font size
lw = 1.5; % line width
mo2d = setup_multi_obj_rand2d(N, M);
mo2d_pareto_points = pareto_opt2d(mo2d);
mo2d_pareto_front   = pareto_front_conts2d(mo2d);
plot_mo_all(1, mo2d, mo2d_pareto_points, mo2d_pareto_front, fs, lw);
end

% ------------------------------------------
% A regular structure constructor to guarantee that
% we have uniform structure names.
function [mo_struct] = make_mo_struct (x1, x2, o1, o2)
  mo_struct = struct('x1', x1, 'x2', x2, 'o1', o1, 'o2', o2);
end

% ----------------------
% setup_multi_obj_rand2d():
%  Sets up all of the objective using random functions
%  In this special structure "multi_obj2d", arrays
%  are created using meshgrid() and stored regularly.
function [mo2d] = setup_multi_obj_rand2d(N, M)
  
% Uniformly sample the parameter space
x1_range = linspace(0, 1, N);
x2_range = linspace(0, 1, M);

% Create a uniform grid of the parameters space based on the range:
[x1, x2] = meshgrid(x1_range, x2_range);
  
% For each one of the parameter values create random values:
% (assignment wants 10 random points, so just make 1-by-N arrays)
f1 = rand(1, N);  % time (to minimize)
f2 = rand(1, N);  % accuracy (to maximize)
x1 = rand(1, N);
x2 = rand(1, N);

 % Put everything in a nice multi-objective optimization structure: 
 mo2d = make_mo_struct(x1, x2, f1, f2);
end

% -------------
% pareto_opt2d()
%   Compute the Pareto-optimal points in 2d.
%   Store the results in a parallel structure of 1-D arrays.
function [mo2d_pareto_points] = pareto_opt2d(mo2d)

% Initialize the results:
x1 = [];
x2 = [];
o1 = [];
o2 = [];

  % Go through a loop all possible points to identify optimal cases    
  [N, M] = size(mo2d.x1);
  for i=1:N
      for j=1:M
          o1_val = mo2d.o1 (i,j);  % time
          o2_val = mo2d.o2 (i,j);  % accuracy
          
          % Compute the number of non-zero elements that are better using
          % nnz(.)
          opt_size = nnz((mo2d.o2 >= o2_val) .* (mo2d.o1 <  o1_val)) + ...
                     nnz((mo2d.o2 >  o2_val) .* (mo2d.o1 <= o1_val));
          
          % Nobody can beat this point. So it is optimal:
          if (opt_size == 0) 
              x1 = [x1, mo2d.x1(i,j)];
              x2 = [x2, mo2d.x2(i,j)];
              o1 = [o1, o1_val];
              o2 = [o2, o2_val];
          end          
      end
  end

  % Create the return structure:
  mo2d_pareto_points = make_mo_struct(x1, x2, o1, o2); 
end

% --------------
% pareto_front2d()
%    Generates a Pareto front by computing the Pareto points
%    and then sorting through them in the right order.
function [mo2d_pareto_front] = pareto_front_conts2d(mo2d)

% Generate the optimal points:
[mo2d_pareto_points] = pareto_opt2d(mo2d);
 
% Sort through the data points:
% Start from the highest value of o1 and the lowest of o2.
% First, re-index based on a descending order from o1.
[o1_sorted, o1_ind] = sort(mo2d_pareto_points.o1, 'ascend');
o2_vals = mo2d_pareto_points.o2 (o1_ind);
x1_vals = mo2d_pareto_points.x1 (o1_ind);
x2_vals = mo2d_pareto_points.x2 (o1_ind);

% Assuming that the objectives are continuous, we do not need to do
% anything else!

% Store the contents in a structure:
mo2d_pareto_front = make_mo_struct(x1_vals, x2_vals, o1_sorted, o2_vals);
end

% --------------
% plot_mo_all():
%   Plots the realization space and Pareto optimal front for visualization.
function plot_mo_all (figure_num, mo_all, mo_points, mo_pareto_front, fs, lw)
  figure(figure_num);
  clf
  hold on
  title('Accuracy vs Time (Pareto front)', 'FontName', 'Arial', 'FontSize', fs);
  xlabel('Time (lower is better)', 'FontName', 'Arial', 'FontSize', fs);
  ylabel('Accuracy (higher is better)', 'FontName', 'Arial', 'FontSize', fs);
  grid on
  axis tight
  plot(mo_all.o1(:), mo_all.o2(:), '.', 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
  h1 = plot(mo_points.o1(:), mo_points.o2(:), 'o', 'MarkerFaceColor', 'b', ...
      'DisplayName', 'Pareto points', 'LineWidth', lw);
  h2 = plot(mo_pareto_front.o1(:), mo_pareto_front.o2(:), 'r-', ...
      'LineWidth', lw, 'DisplayName', 'Pareto front');
  legend([h1 h2], 'Pareto points', 'Pareto front', 'Location', 'northeast', ...
      'FontName', 'Arial', 'FontSize', fs-2);
  set(gca, 'FontName', 'Arial', 'FontSize', fs);
  hold off
  saveas(gcf, 'fig_1.png');
end
