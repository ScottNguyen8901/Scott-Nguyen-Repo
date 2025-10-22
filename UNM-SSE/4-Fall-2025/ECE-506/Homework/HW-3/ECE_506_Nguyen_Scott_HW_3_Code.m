%% Scott Nguyen ECE 506 HW 3
clear;
clc;
close all;

% Common settings
tols_r = [1e-8, 1e-12];
tols_n = [1e-8, 1e-8, 1e-8];
maxit  = 100;

%% -------- Problem 1: Bisection --------
% 1(a)i
bisect_case(@(x) 2*x - 5, 0, 5, 2.5, '1a_i_WorkingExample');
% 1(a)ii
bisect_case(@(x) 2*x + 5, 0, 5, [], '1a_ii_FailureExample');
% 1(b)i
bisect_case(@(x) x.^2 - 5*x + 6, 0, 2.5, 2, '1b_i_Quadratic_Working');
% 1(b)iiA
bisect_case(@(x) x.^2 + 1, -2, 2, [], '1b_ii_Failure_NoSignChange');
% 1(b)iiB
bisect_case(@(x) (x - 2).^2, 1, 3, 2, '1b_ii_Failure_DoubleRoot_NoSignChange');
% 1(c)ii
f_c = @(x) cos(x) - x; x_star = fzero(f_c, [0 1]);
bisect_case(f_c, 0, 1, x_star, '1c_ii_ContinuousFunction');

%% -------- Problem 2: Bisection on f'(x)=0 --------
% 2(a)i
bisect_case(@(x) 2*(x - 2.5), 0, 5, 2.5, '2a_i_fprime_Linear_BracketedMin');
% 2(a)ii
bisect_case(@(x) 2*(x + 5), 0, 5, [], '2a_ii_fprime_Linear_NoSignChange_Fail');
% 2(b)i
bisect_case(@(x) 3*(x.^2 - 4*x + 3), 2.5, 3.5, 3, '2b_i_fprime_Quadratic_BracketedMin');
% 2(b)iiA
bisect_case(@(x) 3*x.^2, -1, 1, [], '2b_iiA_fprime_Quadratic_NoSignChange_Fail');
% 2(b)iiB
bisect_case(@(x) 3*x.^2, -0.5, 0.5, 0, '2b_iiB_fprime_Quadratic_DoubleRoot_NoSignChange');
% 2(c)ii
bisect_case(@(x) -sin(x), 2, 4, pi, '2c_ii_fprime_Continuous_MinExample');

%% -------- Problem 3: Newton --------
% 3(a)i
newton_case(@(x) 2*x - 5, @(x) 2, 0, [0 5], 2.5, '3a_i_Newton_Linear_Working');
% 3(a)ii
newton_case(@(x) 1 + 0*x, @(x) 0*x, 0, [-1 1], [], '3a_ii_Newton_Linear_Fail');
% 3(b)i
newton_case(@(x) x.^2 - 5*x + 6, @(x) 2*x - 5, 0, [-1 3], 2, '3b_i_Newton_Quadratic_Working');
% 3(b)ii
newton_case(@(x) (x - 2).^2 + 1, @(x) 2*(x - 2), 0, [-1 4], [], '3b_ii_Newton_Quadratic_NoRealRoot_Fail');
% 3(c)i
newton_case(@(x) sin(x), @(x) cos(x), 0.2, [-0.5 0.5], 0, '3c_i_Newton_Sin_Working');
% 3(c)ii-B
newton_case(@(x) x.^3 - 2*x + 2, @(x) 3*x.^2 - 2, sqrt(2/3), [-2 2], [], '3c_ii_Newton_DerivZero_Breakdown');
% 3(c)iii (DS on same cubic for real root)
f3 = @(x) x.^3 - 2*x + 2; x_ref3 = fzero(f3, [-3 -1]);
ds_case(f3, 0.0, [-3 1], x_ref3, '3c_iii_ds_Cubic_SameFunction_RealRoot');

%% -------- Problem 4: Newton on f'(x)=0 --------
% 4(a)i
f = @(x) (x - 2.5).^2; df = @(x) 2*(x - 2.5); d2f = @(x) 2 + 0*x;
newton_fprime_case(df, d2f, 0, [0 5], 2.5, '4a_i_NewtonOnGrad_Linear_Working');
% 4(a)ii
f = @(x) x + 1; df = @(x) 1 + 0*x; d2f = @(x) 0*x;
newton_fprime_case(df, d2f, 0, [-5 5], [], '4a_ii_NewtonOnGrad_Linear_Fail');
% 4(b)i
f = @(x) x.^3 - 6*x.^2 + 9*x; df = @(x) 3*(x.^2 - 4*x + 3); d2f = @(x) 6*x - 12;
newton_fprime_case(df, d2f, 3.2, [2.5 3.5], 3, '4b_i_NewtonOnGrad_Quadratic_Working');
% 4(b)ii
f = @(x) x.^3; df = @(x) 3*x.^2; d2f = @(x) 6*x;
newton_fprime_case(df, d2f, 0, [-1 1], [], '4b_ii_NewtonOnGrad_DerivZeroAtStart');
% 4(c)i
f = @(x) cos(x); df = @(x) -sin(x); d2f = @(x) -cos(x);
newton_fprime_case(df, d2f, 3.2, [2 4], pi, '4c_i_NewtonOnGrad_Smooth_Working');
% 4(c)ii
f = @(x) x.^3; df = @(x) 3*x.^2; d2f = @(x) 6*x;
newton_fprime_case(df, d2f, 0, [-0.5 0.5], [], '4c_ii_NewtonOnGrad_FlatStationary_Fail');
% 4(c)iii (DS on gradient)
f = @(x) x.^3 - 2*x + 2; df = @(x) 3*x.^2 - 2; x_ref = fzero(df, [0 1]);
ds_fprime_case(df, 0.0, [-2 2], x_ref, '4c_iii_DS_OnGrad_Working');

%% -------- Problem 5: Comparative Study (f'(x)=0) --------
f  = @(x) (x - 2).^2 + sin(x);
df = @(x) 2*(x - 2) + cos(x);
d2f = @(x) 2 - sin(x);
true_min = fminbnd(f, 0, 4);

res = struct('name',{},'x_hist',{},'fe_hist',{},'x_final',{},'f_final',{}, ...
             'df_final',{},'iters',{},'fe_total',{},'time_sec',{},'pct_err',{});

tic; [x_b, fe_b] = bisection(0, 4, tols_r, df, 200); t_b = toc;
res(end+1) = make_res('Bisection', x_b, fe_b, f, df, true_min, t_b);

tic; [x_n, fe_n] = newton(1.0, 1, tols_n, df, d2f, 200); t_n = toc;
res(end+1) = make_res('Newton', x_n, fe_n, f, df, true_min, t_n);

tic; [x_d, fe_d] = ds_method(1.0, 1, tols_n, df, 200); t_d = toc;
res(end+1) = make_res('DS', x_d, fe_d, f, df, true_min, t_d);

compare_fx_iters_problem5(f, res, [0 4], true_min, '5a_fx_vs_iters');
print_summary_table_problem5(res, true_min);

%% ---------------- Helpers ----------------
function bisect_case(f, a, b, x_star, ttl)
    [x, fe] = bisection(a, b, [1e-8, 1e-12], f, 100);
    visualize_rootfinding_results(f, a, b, x, fe, x_star, ttl);
end

function newton_case(f, df, x0, ab, x_star, ttl)
    [x, fe] = newton(x0, 1, [1e-8, 1e-8, 1e-8], f, df, 100);
    visualize_rootfinding_results(f, ab(1), ab(2), x, fe, x_star, ttl);
end

function ds_case(f, x0, ab, x_star, ttl)
    [x, fe] = ds_method(x0, 1, [1e-8, 1e-8, 1e-8], f, 100);
    visualize_rootfinding_results(f, ab(1), ab(2), x, fe, x_star, ttl);
end

function newton_fprime_case(g, gp, x0, ab, x_star, ttl)
    [x, fe] = newton(x0, 1, [1e-8, 1e-8, 1e-8], g, gp, 100);
    visualize_rootfinding_results(g, ab(1), ab(2), x, fe, x_star, ttl);
end

function ds_fprime_case(g, x0, ab, x_star, ttl)
    [x, fe] = ds_method(x0, 1, [1e-8, 1e-8, 1e-8], g, 100);
    visualize_rootfinding_results(g, ab(1), ab(2), x, fe, x_star, ttl);
end

function [T, h] = visualize_rootfinding_results(f, a, b, x, fe, x_star, ttl)
    if nargin < 6, x_star = []; end
    if nargin < 7, ttl = ''; end
    set(groot, 'DefaultAxesFontSize', 14, 'DefaultLineLineWidth', 2);

    y  = f(x); fa = f(a); fb = f(b); xh = x(end);
    h = figure; hold on; grid on; box on;
    xp = linspace(a, b, 200);
    plot(xp, f(xp), 'k-', 'DisplayName', 'f or f'''); yline(0, '--k', 'LineWidth', 1);
    plot(a, fa, 'go', 'MarkerFaceColor', 'g');
    plot(b, fb, 'ro', 'MarkerFaceColor', 'r');
    plot(x, y, 'bo-', 'MarkerFaceColor', 'b');
    plot(xh, f(xh), 'ms', 'MarkerFaceColor', 'm');
    xlabel('x'); ylabel('value'); title(strrep(ttl, '_', ' '));
    legend('Location','best'); hold off;

    if ~isempty(ttl)
        if ~exist('plots', 'dir'), mkdir('plots'); end
        saveas(h, fullfile('plots', [ttl '.png']));
    end

    it = (0:numel(x)-1).';
    if isempty(x_star)
        T = table(it, x(:), y(:), fe(:), 'VariableNames', {'Iteration','x_k','y_k','FunEvals'});
    else
        T = table(it, x(:), y(:), fe(:), abs(x(:) - x_star), ...
                  'VariableNames', {'Iteration','x_k','y_k','FunEvals','AbsError'});
    end
    disp(' '); disp('Table of Iteration Results:'); disp(T);
end

%% ----- Problem 5 helpers -----
function R = make_res(name, xhist, fehist, f, df, x_true, t)
    x_final = xhist(end);
    R = struct( ...
      'name', name, 'x_hist', xhist(:), 'fe_hist', fehist(:), ...
      'x_final', x_final, 'f_final', f(x_final), 'df_final', df(x_final), ...
      'iters', numel(xhist)-1, ...
      'fe_total', tern(isempty(fehist), NaN, fehist(end)), ...
      'time_sec', t, ...
      'pct_err', 100*abs(x_final - x_true)/max(1,abs(x_true)));
end

function out = tern(cond, a, b)
    if cond, out = a; else, out = b; end
end

function compare_fx_iters_problem5(f, res, ab, true_min, ttl)
    colors = lines(numel(res));
    x_plot = linspace(ab(1), ab(2), 400); y_plot = f(x_plot);
    figure('Name', ttl, 'Position', [100 100 1300 400]);
    for i = 1:numel(res)
        subplot(1,3,i); hold on; grid on; box on;
        plot(x_plot, y_plot, 'k-', 'LineWidth', 1.5, 'DisplayName','f(x)');
        plot(res(i).x_hist, f(res(i).x_hist), 'o-', ...
            'Color', colors(i,:), 'MarkerFaceColor', colors(i,:), ...
            'DisplayName', sprintf('%s iterations', res(i).name));
        xline(true_min, '--', 'Color', [0.5 0 0.5], 'DisplayName','True minimum');
        xlabel('x'); ylabel('f(x)'); title(sprintf('%s Method', res(i).name));
        legend('Location','best');
    end
    if ~isempty(ttl)
        if ~exist('plots','dir'), mkdir('plots'); end
        saveas(gcf, fullfile('plots',[ttl '.png']));
    end
end

function print_summary_table_problem5(res, x_true)
    n = numel(res);
    T = table( ...
      string({res.name})', ...
      [res.x_final]', [res.f_final]', [res.iters]', ...
      [res.fe_total]', [res.time_sec]', [res.pct_err]', ...
      'VariableNames', {'Method','x_final','f_x_final','Iterations','FunEvals','Time_sec','PctErr_from_true_x'});
    disp(' '); disp('Problem 5 Summary:'); disp(T);
end