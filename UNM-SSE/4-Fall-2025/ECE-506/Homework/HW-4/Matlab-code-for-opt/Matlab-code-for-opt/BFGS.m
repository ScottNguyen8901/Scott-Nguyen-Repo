function [x, x_all, stats] = BFGS (fun, G_fun, initial_x, initial_H, vec_eps)
% [x, x_all, stats] = BFGS (fun, G_fun, initial_x, initial_H, vec_eps)
% BFGS (inverse-Hessian form) with Wolfe line search.
%
% Inputs:
%   fun        : f(x) -> scalar
%   G_fun      : g(x) -> n×1
%   initial_x  : n×1
%   initial_H  : n×n (initial inverse-Hessian guess)
%   vec_eps    : gradient norm tolerance
%
% Outputs:
%   x          : final estimate
%   x_all      : path of iterates (n×(k+1))
%   stats      : struct with per-iteration metrics:
%                .f_evals(k), .g_evals(k), .H_evals(k)=0, .mem(k)=n^2+3n

% ---------- Line-search parameters ----------
initial_a_i = 1.0;
a_max = 20;
c1    = 1.0e-4;   % Armijo
c2    = 0.9;      % curvature (strong Wolfe)

% ---------- Initialization ----------
x_k = initial_x(:);
H_k = initial_H;
g_k = G_fun(x_k);

Max_iterations = 1000;

N     = length(x_k);
Ident = eye(N);

x_all = x_k;

% ---------- Instrumentation ----------
C     = EvalCounter();
f_cnt = counted_f(fun,  C);
g_cnt = counted_g(G_fun, C);

stats.f_evals = [];
stats.g_evals = [];
stats.H_evals = [];   % not used in BFGS
stats.mem     = [];

k = 0;
while (norm(g_k) > vec_eps) && (k < Max_iterations)
    % Reset counts at the start of this iteration
    C.reset();

    % Search direction
    p_k = - H_k * g_k;

    % Wolfe line search (counted handles)
    f_k     = f_cnt(x_k);
    phi     = @(a) f_cnt(x_k + a*p_k);
    phi_der = @(a) (g_cnt(x_k + a*p_k).') * p_k;  % (A.16) NW p.583

    a_opt = line_search(phi, phi_der, f_k, (g_k.' * p_k), ...
                        initial_a_i, a_max, c1, c2);

    % Step and store
    x_k_plus_1 = x_k + a_opt * p_k;
    x_all = [x_all, x_k_plus_1];

    % Gradient update (counted)
    g_k_plus_1 = g_cnt(x_k_plus_1);

    % BFGS update
    s_k = x_k_plus_1 - x_k;
    y_k = g_k_plus_1 - g_k;

    sty = (y_k.') * s_k;
    if sty <= 0
        % y_k^T s_k must be positive for BFGS; exit gracefully
        % (see Nocedal & Wright §6.1)
        x_k = x_k_plus_1;
        g_k = g_k_plus_1;
        % log counts/memory before breaking
        n = numel(x_k);
        stats.f_evals(end+1,1) = C.f;
        stats.g_evals(end+1,1) = C.g;
        stats.H_evals(end+1,1) = 0;
        stats.mem(end+1,1)     = n^2 + 3*n;
        break;
    end

    rho_k = 1 / sty;
    A     = Ident - rho_k * (s_k * y_k.');
    A_t   = Ident - rho_k * (y_k * s_k.');
    H_k   = A * H_k * A_t + rho_k * (s_k * s_k.');

    % Advance
    x_k = x_k_plus_1;
    g_k = g_k_plus_1;
    k   = k + 1;

    % Log per-iteration evals and essential memory (in floats)
    n = numel(x_k);
    stats.f_evals(end+1,1) = C.f;
    stats.g_evals(end+1,1) = C.g;
    stats.H_evals(end+1,1) = 0;
    stats.mem(end+1,1)     = n^2 + 3*n;  % x(n) + g(n) + p(n) + H(n^2)
end

x = x_k;
end