function [] = plot_fun_opt (str, f, x1, x2, num_of_lines)
% plot_fun_opt (str, f, x1, x2, num_of_lines)
% 
% Inputs:
%  str:    A title string for the plot.
%  f:      The function to plot.
%  x1, x2: The coordinate system to plot.
%  num_of_lines: number of contour lines.
%
% Outputs:
%  A contour graph of f and a surf plot of f.
%
% Example:
%  N1 = 10;
%  N2 = 10;
%  [x1, x2] = meshgrid(0:N1-1, 0:N2-1);
%  f = @(x) 2*(x(1)-N1/2).^2 + (x(2)-N2/2).^2;
%  plot_fun_opt ('Quad', f, x1, x2, 20);

% Build the function values.
fvals = zeros(size(x1)); 
for i=1:length(x1)
    for j=1:length(x2)
        fvals(i,j) = f([x1(i,j), x2(i,j)]);
    end
end

% Plot them.
figure;
[C, h] = contour(x1, x2, fvals, num_of_lines);
clabel(C,h);
title(str);

figure;
surf(x1, x2, fvals);
title(str);
end

