function [x, x_all, stats] = newton (fun, g_fun, Hessian_fun, initial_x, vec_eps)
% [x, x_all, stats] = newton (fun, g_fun, Hessian_fun, initial_x, vec_eps)
% Newton's method with backtracking line search (Nocedal & Wright Proc. 3.1).
%
% Inputs:
%   fun         : f(x) -> scalar
%   g_fun       : g(x) -> n×1
%   Hessian_fun : H(x) -> n×n
%   initial_x   : n×1 initial guess
%   vec_eps     : gradient-norm tolerance
%
% Outputs:
%   x      : final estimate
%   x_all  : path of iterates (n×(k+1))
%   stats  : struct with per-iteration metrics:
%            .f_evals(k), .g_evals(k), .H_evals(k), .mem(k) ~ n^2 + 3n
%
% Notes:
% - Uses counted wrappers so ALL evaluations (incl. line search) are tracked.
% - Essential memory only (exclude logging/plots).

% ---- Line-search params ----
initial_alpha = 1.0;
c   = 1e-4;
rho = 0.5;

max_iterations = 1000;

% ---- Init ----
x_k   = initial_x(:);
x_all = x_k;

% ---- Instrumentation: counters & counted wrappers ----
C     = EvalCounter();         % requires EvalCounter.m
f_cnt = counted_f(fun, C);     % see helper wrappers counted_f/g/H
g_cnt = counted_g(g_fun, C);
H_cnt = counted_H(Hessian_fun, C);

stats.f_evals = [];
stats.g_evals = [];
stats.H_evals = [];
stats.mem     = [];

k = 0;
while true
    if k >= max_iterations
        break;
    end

    % Reset per-iteration counters
    C.reset();

    % Gradient & stopping test
    g = g_cnt(x_k);
    if norm(g) < vec_eps
        break;
    end

    % Hessian and Newton step (naive solve assumes PD)
    H = H_cnt(x_k);
    warning('off', 'MATLAB:singularMatrix');
    warning('off', 'MATLAB:nearlySingularMatrix');

    p_k = - H \ g;

    % Backtracking line search (use COUNTED f)
    f_x_k   = f_cnt(x_k);
    alpha_k = backtrack_line_search (f_cnt, f_x_k, x_k, p_k, g, initial_alpha, c, rho);

    % Step and store
    x_k = x_k + alpha_k * p_k;
    x_all(:, end+1) = x_k;

    % Log counts & essential memory
    n = numel(x_k);
    stats.f_evals(end+1,1) = C.f;
    stats.g_evals(end+1,1) = C.g;
    stats.H_evals(end+1,1) = C.H;
    stats.mem(end+1,1)     = n^2 + 3*n;  % x(n) + g(n) + p(n) + H(n^2)

    k = k + 1;
end

x = x_k;
end