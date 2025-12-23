%% test_lambert_COMPARE_Izzo_vs_Oldenhuis_with_propagateOrbit.m
clear; clc;

% Use SI units because Aerospace Toolbox functions are in m and m/s
mu = 3.986004418e14;  % Earth mu [m^3/s^2]

%% 1) Define a "truth" orbit using Keplerian elements (degrees) and get epoch state
a_m   = 11000e3;   % [m]
ecc   = 0.1;
incl  = 2;        % [deg]
RAAN  = 2;        % [deg]
argp  = 2;        % [deg]
nu0   = 1;        % [deg]

% keplerian2ijk returns r [m], v [m/s]
[rEpoch, vEpoch] = keplerian2ijk(a_m, ecc, incl, RAAN, argp, nu0);

%% 2) Propagate the truth orbit using Aerospace Toolbox propagateOrbit
tof   = 2.5*3600;                    % [s]
epoch = datetime(2025,1,1,12,0,0);   % arbitrary epoch
dt    = 60;                          % sample cadence [s]
time  = epoch + seconds(0:dt:tof);   % time vector (epoch is time(1))

[posTruth, velTruth] = propagateOrbit( ...
    time, rEpoch, vEpoch, ...
    PropModel="two-body-keplerian");

posTruth = reshape(squeeze(posTruth), 3, []);
velTruth = reshape(squeeze(velTruth), 3, []);

r1       = posTruth(:,1);
v1_truth = velTruth(:,1);
r2       = posTruth(:,end);
v2_truth = velTruth(:,end);

disp("=== Truth Summary ===")
TruthSummary = table( ...
    norm(r1)/1e3, norm(r2)/1e3, tof, dt, ...
    'VariableNames', {'|r1|_km','|r2|_km','TOF_s','dt_s'});
disp(TruthSummary)

%% 3) Solve Lambert using YOUR Izzo solver (SI units)
solS = lambertIzzo(r1, r2, tof, mu, 'longway', false);
solL = lambertIzzo(r1, r2, tof, mu, 'longway', true);

if ~isempty(solS.v1)
    solI = solS; wayI = "short-way";
elseif ~isempty(solL.v1)
    solI = solL; wayI = "long-way";
else
    error("Izzo solver returned no solutions.");
end

% Prefer M=0 if present
idxI = find(solI.M == 0, 1, 'first');
if isempty(idxI), idxI = 1; end

v1_izzo = solI.v1(idxI,:).';
v2_izzo = solI.v2(idxI,:).';

%% 4) Propagate Izzo solution to get r2 and v2 at tf
[posIzzo, velIzzo] = propagateOrbit(time, r1, v1_izzo, PropModel="two-body-keplerian");
posIzzo = reshape(squeeze(posIzzo), 3, []);
velIzzo = reshape(squeeze(velIzzo), 3, []);
r2_izzo      = posIzzo(:,end);
v2_izzo_prop = velIzzo(:,end);

%% 5) Solve Lambert using Oldenhuis lambert.m (FileExchange)
% Oldenhuis expects: r[km], tf[days], mu[km^3/s^2], m integer
mu_km   = mu/1e9;
r1_km   = (r1/1e3).';
r2_km   = (r2/1e3).';
tf_days = tof/86400;
m_rev   = 0;

% Both ways available; match Izzo way
[V1oS_kmps, V2oS_kmps, ~, flagS] = lambert(r1_km, r2_km, +tf_days, +m_rev, mu_km);
[V1oL_kmps, V2oL_kmps, ~, flagL] = lambert(r1_km, r2_km, -tf_days, +m_rev, mu_km);

v1_oldS = V1oS_kmps(:) * 1e3;  v2_oldS = V2oS_kmps(:) * 1e3;
v1_oldL = V1oL_kmps(:) * 1e3;  v2_oldL = V2oL_kmps(:) * 1e3;

if wayI == "short-way"
    v1_old = v1_oldS; v2_old = v2_oldS; flagOld = flagS; wayOld = "short-way";
else
    v1_old = v1_oldL; v2_old = v2_oldL; flagOld = flagL; wayOld = "long-way";
end

if flagOld ~= 1
    warning("Oldenhuis lambert.m returned exitflag=%d (may indicate no solution).", flagOld);
end

%% 6) Propagate Oldenhuis solution to get r2 and v2 at tf
[posOld, velOld] = propagateOrbit(time, r1, v1_old, PropModel="two-body-keplerian");
posOld = reshape(squeeze(posOld), 3, []);
velOld = reshape(squeeze(velOld), 3, []);
r2_old      = posOld(:,end);
v2_old_prop = velOld(:,end);

%% 7) Build ONE comparison table (Izzo vs Oldenhuis)
% Position closure at tf (propagated)
dr2_izzo = r2_izzo - r2;
dr2_old  = r2_old  - r2;

% Velocity mismatch vs truth
dv1_izzo_truth = v1_izzo - v1_truth;
dv2_izzo_truth = v2_izzo_prop - v2_truth;

dv1_old_truth  = v1_old - v1_truth;
dv2_old_truth  = v2_old_prop - v2_truth;

% Solver-to-solver differences
dv1_solver = v1_izzo - v1_old;
dv2_solver = v2_izzo - v2_old;

CompareTable = table( ...
    wayI, solI.M(idxI), solI.branch(idxI), solI.x(idxI), solI.y(idxI), ...                 % Izzo metadata
    wayOld, flagOld, ...                                                                  % Old metadata
    norm(dr2_izzo), norm(dr2_old), ...                                                    % position closure at tf
    norm(dv1_izzo_truth), norm(dv1_old_truth), ...                                        % v1 vs truth
    norm(dv2_izzo_truth), norm(dv2_old_truth), ...                                        % v2 (prop) vs truth
    norm(dv1_solver), norm(dv2_solver), ...                                               % solver-to-solver
    'VariableNames', { ...
        'Izzo_Way','Izzo_M','Izzo_Branch','Izzo_x','Izzo_y', ...
        'Old_Way','Old_exitflag', ...
        '||Δr2||_Izzo_m','||Δr2||_Old_m', ...
        '||Δv1(Izzo-Truth)||_mps','||Δv1(Old-Truth)||_mps', ...
        '||Δv2prop(Izzo-Truth)||_mps','||Δv2prop(Old-Truth)||_mps', ...
        '||Δv1(Izzo-Old)||_mps','||Δv2(Izzo-Old)||_mps'});

disp("=== ONE TABLE: Izzo vs Oldenhuis (with propagated endpoint) ===")
disp(CompareTable)

%% 8) (Optional) component-wise check at tf (position + velocity)
% Uncomment if you want detailed components too.
% Component = ["x";"y";"z"];
% CompareComponents = table( ...
%     Component, ...
%     dr2_izzo(:), dr2_old(:), ...
%     dv2_izzo_truth(:), dv2_old_truth(:), ...
%     'VariableNames', {'Component','dr2_Izzo_m','dr2_Old_m','dv2prop_Izzo-Truth_mps','dv2prop_Old-Truth_mps'});
% disp("=== Component-wise errors at tf ===")
% disp(CompareComponents)
