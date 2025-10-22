% %% ECE 506 - HW4 Problems 1(c)(ii) and 1(e)
% % Author: Scott Nguyen | Fall 2025
% 
% clear; clc; close all;
% 
% % ---------- Global Plot Settings ----------
% set(groot, 'defaultAxesFontName', 'Times New Roman');
% set(groot, 'defaultTextFontName', 'Times New Roman');
% set(groot, 'defaultAxesFontSize', 12);
% set(groot, 'defaultTextFontSize', 12);
% set(groot, 'defaultLineLineWidth', 2);
% 
% % ---------- Ensure output folder ----------
% if ~exist('plots','dir'), mkdir('plots'); end
% 
% %% ---------------- Problem 1(c)(ii): alpha*(a,x1) ----------------
% a_vals  = linspace(-5, 5, 200);
% x1_vals = linspace(-2, 2, 200);
% [A, X1] = meshgrid(a_vals, x1_vals);
% 
% % alpha* = 1/(2a); undefined at a = 0
% Alpha = 1 ./ (2*A);
% Alpha(abs(A) < 1e-8) = NaN;
% 
% figure;
% surf(A, X1, Alpha, 'EdgeColor', 'none');
% xlabel('a'); ylabel('x_1'); zlabel('\alpha^*');
% title('\bfStepsize \alpha^* as a function of (a, x_1)');
% colorbar; view(40,30); grid on;
% saveas(gcf, fullfile('plots','1cii.png'));
% 
% %% ---------------- Problem 1(e): verify 1(c) with provided code ----------------
% % 1D quadratic: f(x) = a x^2, gradient = 2 a x
% a = 2;
% fun   = @(x) a*x^2;       % scalar input x
% g_fun = @(x) 2*a*x;
% 
% x_0     = 1;              % initial guess
% vec_eps = 1e-6;           % stopping threshold
% 
% % Use quiet SD to avoid alpha_k prints
% [x_min, x_all] = steepest_descent_quiet(fun, g_fun, x_0, vec_eps);
% fprintf('Computed x* = %.6f (Analytical x* = 0, Error = %.2e)\n', x_min, abs(x_min));
% 
% % Plot f(x) and iteration points
% x_vals = linspace(-2, 2, 200);
% f_vals = arrayfun(fun, x_vals);
% 
% figure; hold on;
% plot(x_vals, f_vals, 'b-');
% plot(x_all(1,:), arrayfun(fun, x_all(1,:)), 'ro-', 'MarkerFaceColor','r');
% xlabel('x_1'); ylabel('f(x)');
% title('\bfVerification of Steepest Descent for f(x) = a x^2');
% legend('f(x)', 'Iteration points', 'Location', 'northwest');
% grid on;
% saveas(gcf, fullfile('plots','1e.png'));

%% Simple run (fixed I/O shapes): each method x each test function once
clear; clc; close all;
syms x y real

%% ---------- Output folder ----------
plotsDir = 'plots';
if ~exist(plotsDir,'dir'); mkdir(plotsDir); end

%% ---------- Symbolic functions ----------
f_sphere = x^2 + y^2;
f_beale  = (1.5 - x + x*y)^2 + (2.25 - x + x*y^2)^2 + (2.625 - x + x*y^3)^2;
f_gold   = (1 + (x + y + 1)^2 * (19 - 14*x + 3*x^2 - 14*y + 6*x*y + 3*y^2)) * ...
           (30 + (2*x - 3*y)^2 * (18 - 32*x + 12*x^2 + 48*y - 36*x*y + 27*y^2));
f_booth  = (x + 2*y - 7)^2 + (2*x + y - 5)^2;

g_sphere = gradient(f_sphere,[x,y]);  H_sphere = hessian(f_sphere,[x,y]);
g_beale  = gradient(f_beale,[x,y]);  H_beale  = hessian(f_beale,[x,y]);
g_gold   = gradient(f_gold,[x,y]);   H_gold   = hessian(f_gold,[x,y]);
g_booth  = gradient(f_booth,[x,y]);  H_booth  = hessian(f_booth,[x,y]);

%% ---------- Numeric wrappers (2x1 input) ----------
fS = matlabFunction(f_sphere,'Vars',{x,y});
gS = matlabFunction(g_sphere,'Vars',{x,y});
HS = matlabFunction(H_sphere,'Vars',{x,y});
F.sphere.f = @(z) fS(z(1),z(2));
F.sphere.g = @(z) gS(z(1),z(2));
F.sphere.H = @(z) HS(z(1),z(2));

fB = matlabFunction(f_beale,'Vars',{x,y});
gB = matlabFunction(g_beale,'Vars',{x,y});
HB = matlabFunction(H_beale,'Vars',{x,y});
F.beale.f  = @(z) fB(z(1),z(2));
F.beale.g  = @(z) gB(z(1),z(2));
F.beale.H  = @(z) HB(z(1),z(2));

fG = matlabFunction(f_gold,'Vars',{x,y});
gG = matlabFunction(g_gold,'Vars',{x,y});
HG = matlabFunction(H_gold,'Vars',{x,y});
F.gold.f   = @(z) fG(z(1),z(2));
F.gold.g   = @(z) gG(z(1),z(2));
F.gold.H   = @(z) HG(z(1),z(2));

fBo = matlabFunction(f_booth,'Vars',{x,y});
gBo = matlabFunction(g_booth,'Vars',{x,y});
HBo = matlabFunction(H_booth,'Vars',{x,y});
F.booth.f  = @(z) fBo(z(1),z(2));
F.booth.g  = @(z) gBo(z(1),z(2));
F.booth.H  = @(z) HBo(z(1),z(2));

%% ---------- Settings ----------
vec_eps   = 1e-6;
initial_H = eye(2);
rng(42);

%% ---------- Problems and methods ----------
problems = {
  struct('name','Sphere','F',F.sphere, 'xstar',[0;0],   'ftrue',0)
  struct('name','Beale','F',F.beale,   'xstar',[3;0.5], 'ftrue',0)
  struct('name','GoldsteinPrice','F',F.gold, 'xstar',[0;-1], 'ftrue',3)
  struct('name','Booth','F',F.booth,   'xstar',[1;3],   'ftrue',0)
};
methods     = {'Steepest Descent','BFGS','Newton'};
method_keys = {'Steepest_Descent','BFGS','Newton'};
Pnum = numel(problems);
M = numel(methods);

%% ---------- Robustness (100 random inits) + stats for Run 1 ----------
N = 100;
metrics = struct();
stats_run1 = struct();
paths_run1 = struct();

for im = 1:M
    mkey = method_keys{im};
    for ip = 1:Pnum
        pname = matlab.lang.makeValidName(problems{ip}.name);
        metrics.(mkey).(pname).iters = zeros(N,1);
        metrics.(mkey).(pname).e1    = zeros(N,1);
        metrics.(mkey).(pname).e2    = zeros(N,1);
    end
end

wb = waitbar(0,'Running optimization trials...');
counter = 0;

for run = 1:N
    for ip = 1:Pnum
        P = problems{ip};
        pname = matlab.lang.makeValidName(P.name);
        r     = rand; theta = 2*pi*rand;
        x0    = P.xstar + r*[cos(theta);sin(theta)];

        [x_sd,path_sd,st_sd]   = steepest_descent(P.F.f,P.F.g,x0,vec_eps);
        [x_bfgs,path_bfgs,st_bfgs] = BFGS(P.F.f,P.F.g,x0,initial_H,vec_eps);
        [x_nt,path_nt,st_nt]   = newton(P.F.f,P.F.g,P.F.H,x0,vec_eps);

        metrics.Steepest_Descent.(pname).iters(run)=numel(st_sd.f_evals);
        metrics.BFGS.(pname).iters(run)=numel(st_bfgs.f_evals);
        metrics.Newton.(pname).iters(run)=numel(st_nt.f_evals);

        metrics.Steepest_Descent.(pname).e1(run)=abs(x_sd(1)-P.xstar(1));
        metrics.Steepest_Descent.(pname).e2(run)=abs(x_sd(2)-P.xstar(2));
        metrics.BFGS.(pname).e1(run)=abs(x_bfgs(1)-P.xstar(1));
        metrics.BFGS.(pname).e2(run)=abs(x_bfgs(2)-P.xstar(2));
        metrics.Newton.(pname).e1(run)=abs(x_nt(1)-P.xstar(1));
        metrics.Newton.(pname).e2(run)=abs(x_nt(2)-P.xstar(2));

        if run==1
            stats_run1.Steepest_Descent.(pname)=st_sd;
            stats_run1.BFGS.(pname)=st_bfgs;
            stats_run1.Newton.(pname)=st_nt;
            paths_run1.Steepest_Descent.(pname)=normalize_path_shape(path_sd);
            paths_run1.BFGS.(pname)=normalize_path_shape(path_bfgs);
            paths_run1.Newton.(pname)=normalize_path_shape(path_nt);
        end
        counter = counter+1;
        waitbar(counter/(N*Pnum),wb);
    end
end
close(wb);

%% ---------- (1) Robustness plot ----------
figure('Color','w');
tl = tiledlayout(Pnum,3,'TileSpacing','compact','Padding','compact');
title(tl,'Robustness across N=100 random initializations');
meth_list = {'Steepest_Descent','BFGS','Newton'};
meth_disp = {'SD','BFGS','Newton'};
cols = lines(3);
for ip=1:Pnum
    pname=matlab.lang.makeValidName(problems{ip}.name);
    nexttile((ip-1)*3+1); hold on; grid on;
    for m=1:3, histogram(metrics.(meth_list{m}).(pname).iters,'Normalization','pdf','DisplayStyle','stairs','EdgeColor',cols(m,:)); end
    ylabel(problems{ip}.name); if ip==1, title('Iterations'); end
    nexttile((ip-1)*3+2); hold on; grid on;
    for m=1:3, histogram(metrics.(meth_list{m}).(pname).e1,'Normalization','pdf','DisplayStyle','stairs','EdgeColor',cols(m,:)); end
    if ip==1, title('|x_1 - x_1^*|'); end
    nexttile((ip-1)*3+3); hold on; grid on;
    for m=1:3, histogram(metrics.(meth_list{m}).(pname).e2,'Normalization','pdf','DisplayStyle','stairs','EdgeColor',cols(m,:)); end
    if ip==1, title('|x_2 - x_2^*|'); end
end
legend(meth_disp,'Orientation','horizontal','NumColumns',3,'Location','southoutside');
exportgraphics(gcf,fullfile(plotsDir,'robustness_combined.png'),'Resolution',300);

%% ---------- (2) Efficiency (Run 1) ----------
figure('Color','w');
tl=tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
title(tl,'Efficiency per Iteration (Run 1)');
labs={'f evals/iter','g evals/iter','H evals/iter','Memory (floats)'};
flds={'f_evals','g_evals','H_evals','mem'};
for p=1:4
    nexttile(p); hold on; grid on;
    for ip=1:Pnum
        pname=matlab.lang.makeValidName(problems{ip}.name);
        for m=1:3
            st=stats_run1.(meth_list{m}).(pname);
            if ~isfield(st,flds{p}),continue;end
            plot(1:numel(st.(flds{p})),st.(flds{p}),'-o','Color',cols(m,:));
        end
    end
    title(labs{p}); xlabel('iteration k');
end
legend(meth_disp,'Location','southoutside','NumColumns',3);
exportgraphics(gcf,fullfile(plotsDir,'eff_run1_2x2.png'),'Resolution',300);

%% ---------- (3) Accuracy (Run 1) ----------
figure('Color','w');
tl=tiledlayout(Pnum,3,'TileSpacing','compact','Padding','compact');
title(tl,'Rate Diagnostics (Run 1): log(e_k), r_1, r_2');
for ip=1:Pnum
    P=problems{ip}; pname=matlab.lang.makeValidName(P.name);
    for m=1:3
        X=paths_run1.(meth_list{m}).(pname);
        if isempty(X), continue; end
        if size(X,1)~=2, X=X.'; end
        err=vecnorm(X-P.xstar,2,1);
        ek=err(1:end-1); ek1=err(2:end);
        nexttile((ip-1)*3+1); hold on; plot(0:numel(ek)-1,log(max(ek,eps)),'-o','Color',cols(m,:));
        nexttile((ip-1)*3+2); hold on; plot(0:numel(ek)-1,ek1./max(ek,eps),'-o','Color',cols(m,:));
        nexttile((ip-1)*3+3); hold on; plot(0:numel(ek)-1,ek1./max(ek.^2,eps),'-o','Color',cols(m,:));
    end
    if ip==1, title('log(e_k)'); title(nexttile((ip-1)*3+2),'r_1'); title(nexttile((ip-1)*3+3),'r_2'); end
end
legend(meth_disp,'Orientation','horizontal','NumColumns',3,'Location','southoutside');
exportgraphics(gcf,fullfile(plotsDir,'accuracy_run1.png'),'Resolution',300);

%% ---------- Helper ----------
function X=normalize_path_shape(X)
if isempty(X),return;end
if size(X,1)~=2 && size(X,2)==2, X=X.';end
end
