function sol = lambertIzzo(r1, r2, tof, mu, varargin)
%LAMBERTIZZO Lambert solver (Izzo 2014) using x-variable + Householder.
%
%   sol = lambertIzzo(r1, r2, tof, mu)
%   sol = lambertIzzo(..., 'longway', true, 'Mmax', 3, 'tol', 1e-12)
%
% Inputs
%   r1, r2 : 3x1 position vectors [length]
%   tof    : time of flight [time] (must be > 0)
%   mu     : gravitational parameter [length^3/time^2] (must be > 0)
%
% Name-Value options
%   'longway' : true/false (default false)  -> choose transfer angle > pi
%   'Mmax'    : max revolutions to attempt (default floor(T*/pi))
%   'tol'     : Householder convergence tolerance on x (default 1e-13)
%   'maxIter' : max Householder iterations (default 20)
%
% Output struct sol with fields:
%   .v1, .v2  : Nx3 arrays of velocity solutions
%   .x, .y    : Nx1 arrays of Lambert x, y values for each solution
%   .M        : Nx1 revolutions count for each solution
%   .branch   : Nx1 string ("single","short","long")
%   .lambda   : scalar lambda for geometry
%   .Tstar    : nondimensional time of flight
%   .s, .c    : geometry scalars

% -------------------- parse options --------------------
p = inputParser;
p.addParameter('longway', false, @(x)islogical(x) && isscalar(x));
p.addParameter('Mmax', [], @(x)isempty(x) || (isscalar(x) && x>=0 && floor(x)==x));
p.addParameter('tol', 1e-13, @(x)isscalar(x) && x>0);
p.addParameter('maxIter', 20, @(x)isscalar(x) && x>=1);
p.parse(varargin{:});
opt = p.Results;

r1 = r1(:); r2 = r2(:);
if numel(r1)~=3 || numel(r2)~=3, error('r1 and r2 must be 3-vectors.'); end
if ~(tof>0) || ~(mu>0), error('tof and mu must be positive.'); end

% -------------------- geometry --------------------
r1n = norm(r1); r2n = norm(r2);
cvec = r2 - r1;
c = norm(cvec);
s = 0.5*(r1n + r2n + c);

% Transfer angle selection (short vs long)
cosTh = dot(r1,r2)/(r1n*r2n);
cosTh = max(-1,min(1,cosTh));
theta = acos(cosTh);
if opt.longway
    theta = 2*pi - theta;
end

% Lancaster-Blanchard lambda (signed)
lambda2 = 1 - c/s;                 % (s-c)/s
lambda2 = max(0, min(1, lambda2));
lambda = sqrt(lambda2);
if theta > pi
    lambda = -lambda;
end

% Unit vectors for reconstruction (Gooding-style)
ir1 = r1 / r1n;
ir2 = r2 / r2n;

ih = cross(ir1, ir2);
ihn = norm(ih);
if ihn < 1e-15
    ih = [0;0;1]; ihn = 1;
end
ih = ih / ihn;

% For longway, flip plane normal so tangential directions are consistent
if opt.longway
    ih = -ih;
end
it1 = cross(ih, ir1);
it2 = cross(ih, ir2);

% Nondimensional time of flight (Izzo)
Tstar = sqrt(2*mu/s^3) * tof;

% Default Mmax guess
if isempty(opt.Mmax)
    Mmax = floor(Tstar/pi);
else
    Mmax = opt.Mmax;
end

% -------------------- solve for all branches --------------------
v1_list = [];
v2_list = [];
x_list  = [];
y_list  = [];
M_list  = [];
branch_list = strings(0,1);

% --- Single revolution (M=0) ---
[x0, ok0] = initial_guess_M0(lambda, Tstar);
if ok0
    [xsol, ysol, ok] = solve_for_x(lambda, Tstar, 0, x0, opt.tol, opt.maxIter);
    if ok
        [v1, v2] = reconstruct_vel(mu, s, c, r1n, r2n, lambda, xsol, ysol, ir1, ir2, it1, it2);
        v1_list = [v1_list; v1.']; %#ok<AGROW>
        v2_list = [v2_list; v2.']; %#ok<AGROW>
        x_list  = [x_list;  xsol]; %#ok<AGROW>
        y_list  = [y_list;  ysol]; %#ok<AGROW>
        M_list  = [M_list;  0];    %#ok<AGROW>
        branch_list(end+1,1) = "single"; %#ok<AGROW>
    end
end

% --- Multi-revolution (M=1..Mmax) ---
for M = 1:Mmax
    [x0l, x0r] = initial_guess_Mgt0(Tstar, M);

    [xL, yL, okL] = solve_for_x(lambda, Tstar, M, x0l, opt.tol, opt.maxIter);
    if okL
        [v1, v2] = reconstruct_vel(mu, s, c, r1n, r2n, lambda, xL, yL, ir1, ir2, it1, it2);
        v1_list = [v1_list; v1.']; %#ok<AGROW>
        v2_list = [v2_list; v2.']; %#ok<AGROW>
        x_list  = [x_list;  xL];   %#ok<AGROW>
        y_list  = [y_list;  yL];   %#ok<AGROW>
        M_list  = [M_list;  M];    %#ok<AGROW>
        branch_list(end+1,1) = "long"; %#ok<AGROW>
    end

    [xS, yS, okS] = solve_for_x(lambda, Tstar, M, x0r, opt.tol, opt.maxIter);
    if okS
        [v1, v2] = reconstruct_vel(mu, s, c, r1n, r2n, lambda, xS, yS, ir1, ir2, it1, it2);
        v1_list = [v1_list; v1.']; %#ok<AGROW>
        v2_list = [v2_list; v2.']; %#ok<AGROW>
        x_list  = [x_list;  xS];   %#ok<AGROW>
        y_list  = [y_list;  yS];   %#ok<AGROW>
        M_list  = [M_list;  M];    %#ok<AGROW>
        branch_list(end+1,1) = "short"; %#ok<AGROW>
    end
end

% -------------------- pack output --------------------
sol = struct();
sol.v1 = v1_list;
sol.v2 = v2_list;
sol.x  = x_list;
sol.y  = y_list;
sol.M  = M_list;
sol.branch = branch_list;
sol.lambda = lambda;
sol.Tstar  = Tstar;
sol.s = s;
sol.c = c;

end

% ======================================================================
% =========================== helpers ==================================
% ======================================================================

function [x0, ok] = initial_guess_M0(lambda, T)
% Izzo Eq. (30) piecewise starter for M=0.
ok = true;
T0 = acos(lambda) + lambda*sqrt(max(0,1-lambda^2));  % Eq. 19 (M=0 at x=0)
T1 = (2/3)*(1 - lambda^3);                           % Eq. 21 (parabolic)
if T <= 0
    ok = false; x0 = NaN; return;
end

if T >= T0
    x0 = (T0/T)^(2/3) - 1;
elseif T < T1
    denom = T*(1 - lambda^5);
    if abs(denom) < 1e-16
        x0 = 2*(T1/T) - 1;
    else
        x0 = (5/2) * (T1*(T1 - T))/denom + 1;
    end
else
    x0 = (T0/T)^(log(2)/log(T1/T0)) - 1;
end

% keep x in admissible-ish range for M=0: [-1, inf)
x0 = max(-0.999999999999, x0);
end

function [x0l, x0r] = initial_guess_Mgt0(T, M)
% Izzo Eq. (31) starters for M>0.
a = ((M*pi + pi) / (8*T))^(2/3);
b = ((8*T) / (M*pi))^(2/3);
x0l = (a - 1) / (a + 1);
x0r = (b - 1) / (b + 1);

% clamp strictly inside (-1,1)
epsx = 1e-12;
x0l = max(-1+epsx, min(1-epsx, x0l));
x0r = max(-1+epsx, min(1-epsx, x0r));
end

function [x, y, ok] = solve_for_x(lambda, Tstar, M, x0, tol, maxIter)
% Householder iterations on f(x)=T(x)-Tstar
x = x0;
ok = false;

for k = 1:maxIter
    [T, d1, d2, d3, y, valid] = tof_and_derivs(lambda, x, M);
    if ~valid || ~isfinite(T)
        return;
    end

    f    = T - Tstar;
    fp   = d1;
    fpp  = d2;
    fppp = d3;

    % Householder (3rd order) update
    A = fp*fp - 0.5*f*fpp;
    B = fp*(fp*fp - f*fpp) + (fppp * f*f)/6;

    if abs(B) < 1e-18 || ~isfinite(B)
        return;
    end

    dx = f*A/B;
    xnew = x - dx;

    if ~isfinite(xnew)
        return;
    end

    if M > 0
        epsx = 1e-12;
        xnew = max(-1+epsx, min(1-epsx, xnew));
    else
        xnew = max(-0.999999999999, xnew);
    end

    if abs(xnew - x) < tol
        x = xnew;
        y = sqrt(1 - lambda^2*(1 - x^2));
        ok = true;
        return;
    end

    x = xnew;
end
end

function [T, dTdx, d2Tdx2, d3Tdx3, y, valid] = tof_and_derivs(lambda, x, M)
% Evaluate nondimensional time of flight T and derivatives (Izzo Eq. 18 & 22),
% with Battin-series fallback near x~1 for M=0.

valid = true;

y2 = 1 - lambda^2*(1 - x^2);
if y2 < 0
    valid = false;
    T = NaN; dTdx = NaN; d2Tdx2 = NaN; d3Tdx3 = NaN; y = NaN; return;
end
y = sqrt(y2);

% Compute T
if M==0 && abs(1 - x) < 1e-3
    T = tof_battin_M0(lambda, x, y);
else
    T = tof_eq18_correct(lambda, x, y, M);
end

den = (1 - x^2);
if abs(den) < 1e-16
    valid = false;
    dTdx = NaN; d2Tdx2 = NaN; d3Tdx3 = NaN;
    return;
end

% Izzo Eq. (22)
dTdx    = (3*T*x - 2 + 2*(lambda^3)*x/y) / den;
d2Tdx2  = (3*T + 5*x*dTdx + 2*(1-lambda^2)*(lambda^3)/(y^3)) / den;
d3Tdx3  = (7*x*d2Tdx2 + 8*dTdx - 6*(1-lambda^2)*(lambda^5)*x/(y^5)) / den;

end

function T = tof_eq18_correct(lambda, x, y, M)
% Correct Izzo Eq. (18):
%   T = (psi + M*pi)/|1-x^2|^(3/2) - (x - lambda*y)/(1-x^2)
%
% psi computed from Eq. (17).

one_minus_x2 = 1 - x^2;

% psi
if x < 1
    srt = sqrt(max(0, one_minus_x2));
    cospsi = x*y + lambda*(one_minus_x2);
    sinpsi = (y - x*lambda)*srt;
    psi = atan2(sinpsi, cospsi);      % (-pi,pi]
    if psi < 0
        psi = psi + 2*pi;
    end
    % fold to [0,pi]
    if psi > pi
        psi = 2*pi - psi;
    end
else
    % hyperbolic (typically only M=0 is meaningful)
    coshpsi = x*y - lambda*(x^2 - 1);
    coshpsi = max(1, coshpsi);
    psi = acosh(coshpsi);
end

abs32 = abs(one_minus_x2)^(3/2);
T = (psi + M*pi)/abs32 - (x - lambda*y)/one_minus_x2;
end

function T = tof_battin_M0(lambda, x, y)
% Battin-series expression used near x~1 for M=0
eta = y - lambda*x;
S1  = 0.5*(1 - lambda - x*eta);

term = 1.0;
sumv = 1.0;
Nmax = 30;
a = 3.0;
c = 2.5; % 5/2
for k = 1:Nmax
    term = term * ((a + (k-1)) / (c + (k-1))) * S1;
    sumv = sumv + term;
    if abs(term) < 1e-16
        break;
    end
end

Q = (4/3) * sumv;
T = 0.5*(eta^3 * Q + 4*lambda*eta);
end

function [v1, v2] = reconstruct_vel(mu, s, c, r1n, r2n, lambda, x, y, ir1, ir2, it1, it2)
% Gooding-style algebraic reconstruction (Izzo Algorithm 1)
gamma = sqrt(mu*s/2);
rho = (r1n - r2n)/c;
sigma = sqrt(max(0, 1 - rho^2));

Vr1 = gamma * ((lambda*y - x) - rho*(lambda*y + x)) / r1n;
Vr2 = -gamma * ((lambda*y - x) + rho*(lambda*y + x)) / r2n;

Vt1 = gamma * sigma * (y + lambda*x) / r1n;
Vt2 = gamma * sigma * (y + lambda*x) / r2n;

v1 = Vr1*ir1 + Vt1*it1;
v2 = Vr2*ir2 + Vt2*it2;
end