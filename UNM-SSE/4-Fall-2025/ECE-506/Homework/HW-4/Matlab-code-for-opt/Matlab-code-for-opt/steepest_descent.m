function [x, x_all, stats] = steepest_descent(fun, g_fun, initial_x, vec_eps)
% [x, x_all, stats] = steepest_descent (fun, g_fun, initial_x, vec_eps)
% Steepest descent with backtracking line-search (Nocedal & Wright, Proc 3.1).
%
% Inputs:
%   fun       : function handle f(x) -> scalar
%   g_fun     : gradient handle g(x) -> n×1
%   initial_x : n×1 initial guess
%   vec_eps   : gradient norm tolerance
%
% Outputs:
%   x         : final estimate
%   x_all     : path (n×(k+1))
%   stats     : struct with per-iteration counts & memory (essential only)
%               .f_evals(k), .g_evals(k), .H_evals(k)=0, .mem(k)=3n

% ---- Line-search parameters ----
initial_alpha = 1.0;
c   = 1e-4;
rho = 0.5;

max_iterations = 1000;

% ---- Init ----
x_k   = initial_x(:);
x_all = x_k;

% ---- Instrumentation: counters & counted wrappers ----
C    = EvalCounter();               % requires EvalCounter.m (handle class)
f_cnt = counted_f(fun, C);          % counted wrappers (see helper funcs)
g_cnt = counted_g(g_fun, C);

stats.f_evals = [];
stats.g_evals = [];
stats.H_evals = [];                 % always zero for SD
stats.mem     = [];

k = 0;
while true
    % Stop if iteration cap hit (protect against pathological line-search loops)
    if k >= max_iterations
        break;
    end

    % Reset per-iteration counters
    C.reset();

    % Gradient & stopping check
    g = g_cnt(x_k);
    if norm(g) < vec_eps
        break;
    end

    % Search direction and line search (use counted f)
    p_k   = -g;
    f_x_k = f_cnt(x_k);
    alpha_k = backtrack_line_search(f_cnt, f_x_k, x_k, p_k, g, initial_alpha, c, rho);

    % Step
    x_k = x_k + alpha_k * p_k;
    x_all(:, end+1) = x_k;

    % Log per-iteration evals and essential memory (in floats)
    n = numel(x_k);
    stats.f_evals(end+1,1) = C.f;
    stats.g_evals(end+1,1) = C.g;
    stats.H_evals(end+1,1) = 0;
    stats.mem(end+1,1)     = 3*n;   % x(n) + g(n) + p(n)

    % Next iter
    k = k + 1;
end

x = x_k;
end
