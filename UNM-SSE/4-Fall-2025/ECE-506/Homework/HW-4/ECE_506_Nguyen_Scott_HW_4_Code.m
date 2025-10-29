%% ECE 506 - HW4 Problems
% Author: Scott Nguyen | Fall 2025

clear
clc
close all

% Global plot settings
set(groot, 'defaultAxesFontName', 'Times New Roman');
set(groot, 'defaultTextFontName', 'Times New Roman');
set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultLineLineWidth', 2);

% Ensure output folder
if ~exist('plots','dir')
    mkdir('plots');
end

%% Problem 1(c)(ii): alpha*(a,x1)
a_vals  = linspace(-5, 5, 200);
x1_vals = linspace(-2, 2, 200);
[A, X1] = meshgrid(a_vals, x1_vals);
Alpha = 1 ./ (2*A);
Alpha(abs(A) < 1e-8) = NaN;

figure;
surf(A, X1, Alpha, 'EdgeColor', 'none');
xlabel('a'); ylabel('x_1'); zlabel('\alpha^*');
title('\bfStepsize \alpha^* as a function of (a, x_1)');
colorbar; view(40,30); grid on;
saveas(gcf, fullfile('plots','1cii.png'));

%% Problem 1(e): verify 1(c)
a = 2;
fun   = @(x) a*x^2;
g_fun = @(x) 2*a*x;
x_0     = 1;
vec_eps = 1e-6;

[x_min, x_all, ~] = steepest_descent(fun, g_fun, x_0, vec_eps);
fprintf('Computed x* = %.6f (Analytical x* = 0, Error = %.2e)\n', x_min, abs(x_min));

x_vals = linspace(-2, 2, 200);
f_vals = arrayfun(fun, x_vals);

figure; hold on;
plot(x_vals, f_vals, 'b-');
plot(x_all(1,:), arrayfun(fun, x_all(1,:)), 'ro-', 'MarkerFaceColor','r');
xlabel('x_1'); ylabel('f(x)');
title('\bfVerification of Steepest Descent for f(x) = a x^2');
legend('f(x)', 'Iteration points', 'Location', 'northwest');
grid on;
saveas(gcf, fullfile('plots','1e.png'));

%% Multi-start benchmarks
clearvars -except Alpha A X1
close all
rng(42)

vec_eps    = 1e-6;
num_starts = 3;
init_box   = [-4 4; -4 4];
n          = 2;
H0         = eye(n);

F.sphere.fun  = @sphere_fun;
F.sphere.grad = @sphere_grad;
F.sphere.hess = @sphere_hess;
Xs.sphere     = [0; 0];

F.beale.fun   = @beale_fun;
F.beale.grad  = @beale_grad;
F.beale.hess  = @(x) num_hess(@beale_grad, x);
Xs.beale      = [3; 0.5];

F.gold.fun    = @gold_fun;
F.gold.grad   = @gold_grad;
F.gold.hess   = @(x) num_hess(@gold_grad, x);
Xs.gold       = [0; -1];

F.booth.fun   = @booth_fun;
F.booth.grad  = @booth_grad;
F.booth.hess  = @booth_hess;
Xs.booth      = [1; 3];

problems = {
  'Sphere',          F.sphere.fun, F.sphere.grad, F.sphere.hess, Xs.sphere;
  'Beale',           F.beale.fun,  F.beale.grad,  F.beale.hess,  Xs.beale ;
  'Goldstein-Price', F.gold.fun,   F.gold.grad,   F.gold.hess,   Xs.gold  ;
  'Booth',           F.booth.fun,  F.booth.grad,  F.booth.hess,  Xs.booth ;
};

% Steepest Descent
Rows = {}; mem_str = '3n';
for i = 1:size(problems,1)
    name  = problems{i,1}; funH  = problems{i,2}; gradH = problems{i,3}; xstar = problems{i,5};
    for s = 1:num_starts
        x0 = [rand_range(init_box); rand_range(init_box)];
        [xf, Xall, counts] = steepest_descent(funH, gradH, x0, vec_eps);
        iters   = size(Xall, 2) - 1;
        f_final = funH(xf);
        gnorm   = norm(gradH(xf));
        conv    = (gnorm <= vec_eps);
        S = rate_stats(Xall, xstar, gradH);
        Rows(end+1,:) = {name, s, x0(:).', xf(:).', f_final, gnorm, iters, ...
                         counts.f, counts.g, counts.H, mem_str, ...
                         S.Rlin_x, S.Rquad_x, char(S.Order_x), S.Rlin_g, S.Rquad_g, char(S.Order_g), conv}; %#ok<AGROW>
    end
end
T_SD = cell2table(Rows, 'VariableNames', table_vars());
disp('=== Steepest Descent ==='); disp(T_SD);

% BFGS
Rows = {}; mem_str = 'n^2 + 3n';
for i = 1:size(problems,1)
    name  = problems{i,1}; funH  = problems{i,2}; gradH = problems{i,3}; xstar = problems{i,5};
    for s = 1:num_starts
        x0 = [rand_range(init_box); rand_range(init_box)];
        [xf, Xall, counts] = BFGS(funH, gradH, x0, H0, vec_eps);
        iters   = size(Xall, 2) - 1;
        f_final = funH(xf);
        gnorm   = norm(gradH(xf));
        conv    = (gnorm <= vec_eps);
        S = rate_stats(Xall, xstar, gradH);
        Rows(end+1,:) = {name, s, x0(:).', xf(:).', f_final, gnorm, iters, ...
                         counts.f, counts.g, counts.H, mem_str, ...
                         S.Rlin_x, S.Rquad_x, char(S.Order_x), S.Rlin_g, S.Rquad_g, char(S.Order_g), conv}; %#ok<AGROW>
    end
end
T_BFGS = cell2table(Rows, 'VariableNames', table_vars());
disp('=== BFGS ==='); disp(T_BFGS);

% Newton
Rows = {}; mem_str = 'n^2 + 2n';
for i = 1:size(problems,1)
    name = problems{i,1}; funH = problems{i,2}; gradH = problems{i,3}; hessH = problems{i,4}; xstar = problems{i,5};
    for s = 1:num_starts
        x0 = [rand_range(init_box); rand_range(init_box)];
        [xf, Xall, counts] = newton(funH, gradH, hessH, x0, vec_eps);
        iters   = size(Xall, 2) - 1;
        f_final = funH(xf);
        gnorm   = norm(gradH(xf));
        conv    = (gnorm <= vec_eps);
        S = rate_stats(Xall, xstar, gradH);
        Rows(end+1,:) = {name, s, x0(:).', xf(:).', f_final, gnorm, iters, ...
                         counts.f, counts.g, counts.H, mem_str, ...
                         S.Rlin_x, S.Rquad_x, char(S.Order_x), S.Rlin_g, S.Rquad_g, char(S.Order_g), conv}; %#ok<AGROW>
    end
end
T_NEW = cell2table(Rows, 'VariableNames', table_vars());
disp('=== Newton ==='); disp(T_NEW);

% Trust-Region (CG-Steihaug)
Rows = {}; mem_str = 'n^2 + 2n';
TR = struct('Delta0',1, 'DeltaMax',1000, 'eta1',0.1, 'eta2',0.75, ...
            'gamma1',0.5, 'gamma2',2.0, 'maxit',200, 'cg_eps',1e-6);
for i = 1:size(problems,1)
    name = problems{i,1}; funH = problems{i,2}; gradH = problems{i,3}; hessH = problems{i,4}; xstar = problems{i,5};
    for s = 1:num_starts
        x0 = [rand_range(init_box); rand_range(init_box)];
        [xf, Xall, Dhist, counts] = tr_cg_driver(funH, gradH, hessH, x0, vec_eps, TR); %#ok<ASGLU>
        iters   = size(Xall, 2) - 1;
        f_final = funH(xf);
        gnorm   = norm(gradH(xf));
        conv    = (gnorm <= vec_eps);
        S = rate_stats(Xall, xstar, gradH);
        Rows(end+1,:) = {name, s, x0(:).', xf(:).', f_final, gnorm, iters, ...
                         counts.f, counts.g, counts.H, mem_str, ...
                         S.Rlin_x, S.Rquad_x, char(S.Order_x), S.Rlin_g, S.Rquad_g, char(S.Order_g), conv}; %#ok<AGROW>
    end
end
T_TRCG = cell2table(Rows, 'VariableNames', table_vars());
disp('=== Trust-Region (CG-Steihaug) ==='); disp(T_TRCG);

% ----------------------- helper functions -----------------------
function v = rand_range(rng2d)
    v = rng2d(1,1) + (rng2d(1,2) - rng2d(1,1)) * rand;
end

function vars = table_vars()
vars = {'Function','RunID','x0','x_final','f_final','grad_norm','Iterations', ...
        'FuncEvals','GradEvals','HessEvals','MemoryPerIter', ...
        'Rlin_error','Rquad_error','Order_error', ...
        'Rlin_grad','Rquad_grad','Order_grad','Converged'};
end

function S = rate_stats(Xall, xs, gradH)
    E = vecnorm(Xall - xs, 2, 1);
    G = arrayfun(@(k) norm(gradH(Xall(:,k))), 1:size(Xall,2));
    rlin_x  = E(2:end) ./ max(E(1:end-1), eps);
    rquad_x = E(2:end) ./ max(E(1:end-1).^2, eps);
    rlin_g  = G(2:end) ./ max(G(1:end-1), eps);
    rquad_g = G(2:end) ./ max(G(1:end-1).^2, eps);
    tail_idx = max(1, ceil(2*numel(rlin_x)/3)) : numel(rlin_x);
    med = @(v) median(v, 'omitnan');
    S.Rlin_x  = med(rlin_x(tail_idx));
    S.Rquad_x = med(rquad_x(tail_idx));
    S.Rlin_g  = med(rlin_g(tail_idx));
    S.Rquad_g = med(rquad_g(tail_idx));
    S.Order_x = classify_order(S.Rlin_x, S.Rquad_x);
    S.Order_g = classify_order(S.Rlin_g, S.Rquad_g);
end

function label = classify_order(Rlin, Rquad)
    if isnan(Rlin)
        label = "none";
        return; 
    end
    if (Rlin < 1e-2) && isfinite(Rquad) && (Rquad > 0) && (Rquad < 1e2)
        label = "quadratic";
        return;
    end
    if (Rlin < 1e-1)
        label = "superlinear";
        return;
    end
    if (Rlin > 0) && (Rlin < 1)
        label = "linear";
    else
        label = "none";
    end
end

function f = sphere_fun(x)
    f = x(1).^2 + x(2).^2;
end

function g = sphere_grad(x)
    g = [2*x(1); 2*x(2)];
end

function H = sphere_hess(~)
    H = [2 0; 0 2];
end

function f = beale_fun(x)
f = (1.5 - x(1) + x(1)*x(2)).^2 + ...
    (2.25 - x(1) + x(1)*x(2).^2).^2 + ...
    (2.625 - x(1) + x(1)*x(2).^3).^2;
end

function g = beale_grad(x)
    A  = 1.5   - x(1) + x(1)*x(2);
    B  = 2.25  - x(1) + x(1)*x(2).^2;
    C  = 2.625 - x(1) + x(1)*x(2).^3;
    dA_dx = -1 + x(2);        dA_dy = x(1);
    dB_dx = -1 + x(2).^2;     dB_dy = 2*x(1)*x(2);
    dC_dx = -1 + x(2).^3;     dC_dy = 3*x(1)*x(2).^2;
    gx = 2*(A*dA_dx + B*dB_dx + C*dC_dx);
    gy = 2*(A*dA_dy + B*dB_dy + C*dC_dy);
    g  = [gx; gy];
end

function f = gold_fun(x)
    u  = x(1) + x(2) + 1;
    v  = 2*x(1) - 3*x(2);
    g1 = 19 - 14*x(1) + 3*x(1).^2 - 14*x(2) + 6*x(1)*x(2) + 3*x(2).^2;
    g2 = 18 - 32*x(1) + 12*x(1).^2 + 48*x(2) - 36*x(1)*x(2) + 27*x(2).^2;
    A  = 1 + u.^2 .* g1;
    B  = 30 + v.^2 .* g2;
    f  = A .* B;
end

function g = gold_grad(x)
    u  = x(1) + x(2) + 1;
    v  = 2*x(1) - 3*x(2);
    g1 = 19 - 14*x(1) + 3*x(1).^2 - 14*x(2) + 6*x(1)*x(2) + 3*x(2).^2;
    g2 = 18 - 32*x(1) + 12*x(1).^2 + 48*x(2) - 36*x(1)*x(2) + 27*x(2).^2;
    dg1_dx = -14 + 6*x(1) + 6*x(2);
    dg1_dy = -14 + 6*x(1) + 6*x(2);
    dg2_dx = -32 + 24*x(1) - 36*x(2);
    dg2_dy = 48 - 36*x(1) + 54*x(2);
    A  = 1 + u.^2 .* g1;
    B  = 30 + v.^2 .* g2;
    dA_dx = 2*u*g1 + u.^2*dg1_dx;
    dA_dy = 2*u*g1 + u.^2*dg1_dy;
    dB_dx = 4*v*g2 + v.^2*dg2_dx;
    dB_dy = -6*v*g2 + v.^2*dg2_dy;
    gx = B .* dA_dx + A .* dB_dx;
    gy = B .* dA_dy + A .* dB_dy;
    g  = [gx; gy];
end

function f = booth_fun(x)
f = (x(1) + 2*x(2) - 7).^2 + (2*x(1) + x(2) - 5).^2;
end

function g = booth_grad(x)
    g1 = 2*(x(1) + 2*x(2) - 7) + 4*(2*x(1) + x(2) - 5);
    g2 = 4*(x(1) + 2*x(2) - 7) + 2*(2*x(1) + x(2) - 5);
    g  = [g1; g2];
end

function H = booth_hess(~)
    H = [10 8; 8 10];
end

function H = num_hess(gradfun, x)
    h   = 1e-6;
    H   = zeros(2,2);
    e1  = [1; 0];
    e2  = [0; 1];
    g1p = gradfun(x + h*e1);
    g1m = gradfun(x - h*e1);
    g2p = gradfun(x + h*e2);
    g2m = gradfun(x - h*e2);
    H(:,1) = (g1p - g1m) / (2*h);
    H(:,2) = (g2p - g2m) / (2*h);
    H = (H + H.') / 2;
end

function [xk, Xall, Dhist, counts] = tr_cg_driver(fun, grad, hess, x0, tol, TR)
    xk    = x0(:);
    Xall  = xk;
    D     = TR.Delta0;
    Dhist = D;
    counts = struct('f',0,'g',0,'H',0);
    for k = 1:TR.maxit
        gk = grad(xk); counts.g = counts.g + 1;
        if norm(gk) <= tol, break; end
        Bk = hess(xk); counts.H = counts.H + 1;
        [p, ~]  = CG_Steihaug(TR.cg_eps, gk, Bk, D);
        mk0 = 0; mkp = gk.'*p + 0.5*p.'*Bk*p;
        denom = -(mkp - mk0); if denom <= 0, denom = eps; end
        fx = fun(xk); counts.f = counts.f + 1;
        fnew = fun(xk + p); counts.f = counts.f + 1;
        rho = (fx - fnew) / denom;
        accepted = false;
        if rho >= TR.eta1
            xk = xk + p; Xall(:, end+1) = xk; %#ok<AGROW>
            accepted = true;
        else
            D = TR.gamma1 * D;
        end
        if rho < 0.25
            D = max(TR.gamma1 * D, 1e-12);
        elseif (rho > TR.eta2) && (abs(norm(p) - D) <= 1e-12)
            D = min(TR.gamma2 * D, TR.DeltaMax);
        end
        Dhist(end+1) = D; %#ok<AGROW>
        if accepted
            gk2 = grad(xk); counts.g = counts.g + 1;
            if norm(gk2) <= tol, break; end
        end
    end
end