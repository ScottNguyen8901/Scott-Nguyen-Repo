function [] = plot_trajectory (name, vec_fun, x_range, y_range, x_all)
%
% plot_trajectory (name, vec_fun, x_range, y_range, x_all)
% Plots the trajectory of iterates in x_all over a 2D
% function given by vec_fun, defined over x_range, y_range.
%
% For an example, see min_rosen.m.
%

[x1, x2] = meshgrid(x_range, y_range);
[cs,h]=contour(x1, x2, vec_fun(x1, x2), 100); hold on;
%clabel(cs,h);
plot(x_all(1,:), x_all(2,:));
plot(x_all(1,:), x_all(2,:), '.');
title(name);
