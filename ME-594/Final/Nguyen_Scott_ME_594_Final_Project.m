clear
close all
clc

%% ------------------------Problem 1---------------------------------------

% Populating constants
constants

% Set parameters
fs = 24;  % Font size
lineWidth = 1.5;  % Line width

% Load the data from the CSV file
data = readmatrix('ME594 centroid data for Final.xlsx');

% Extracting coordinates for each image
x_cam = data(:, [1, 4, 7, 10, 13]);  % x-coordinates for images 1-5
y_cam = data(:, [2, 5, 8, 11, 14]);  % y-coordinates for images 1-5
colors = ['r', 'b', 'g', 'm', 'c'];  % Colors for the images

% Plotting the coordinates for each camera set
figure;
hold on;
for i = 1:5
    plot(x_cam(:,i), y_cam(:,i), [colors(i), 'o'], 'DisplayName', sprintf('Image %d', i), 'LineWidth', lineWidth);
end

% Add title, labels, and legend
title('Camera Coordinates for Images 1-5', 'FontSize', fs, 'FontWeight', 'bold');
xlabel('X Coordinates (pixels)', 'FontSize', fs, 'FontWeight', 'bold');
ylabel('Y Coordinates (pixels)', 'FontSize', fs, 'FontWeight', 'bold');
legend show;

% Enable grid and set axis properties
grid on;
set(gca, 'FontSize', fs, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

hold off;

%% ------------------------Problem 2---------------------------------------

rng(464876754);
% rng('shuffle');

% Initial setup
date0 = datetime('2025-05-05 03:18:34');
L = deg2rad(35.0844);        % Albuquerque latitude
lon = deg2rad(-106.6504);    % Albuquerque longitude

% Initial ISS state in ECI [km; km/s]
r0 = [5165.8127; -1938.2678; -3951.6672];
v0 = [4.7525; 4.4652; 4.0249];
x_ISS_ECI_0 = [r0; v0];

% Time span and integration
tspan = 0:60*21;
opts = odeset('RelTol',1e-12, 'AbsTol',1e-12);
[t, x_ISS_ECI] = ode45(@(t,y) two_body_ode(t, y), tspan, x_ISS_ECI_0, opts);
r_ISS_ECI = x_ISS_ECI(:,1:3);

% Allocate arrays
N = length(t);
RA_Dec_true = zeros(N, 2); 
RA_Dec_pert = zeros(N, 2);
Az_El_true  = zeros(N, 2); 
Az_El_pert  = zeros(N, 2);

date_vec = date0 + seconds(t);
sigma    = deg2rad(0.01);         % Noise std [rad]

% Main loop
for i = 1:N
    theta_g = GMST(date_vec(i));
    theta = lon + theta_g;
    Q = rot_ijk_sez(L, theta);

    r_th = Q * r_ISS_ECI(i,:)';
    u_th = (r_th - r_s_TH) / norm(r_th - r_s_TH);

    % Azimuth and Elevation
    az = mod(atan2(u_th(1), u_th(2)), 2*pi);
    el = asin(-u_th(3));
    Az_El_true(i,:) = [az, el];
    Az_El_pert(i,:) = [az, el] + sigma * randn(1,2);

    % RA and Dec
    r_s_eci = Q' * r_s_TH;
    u_eci = (r_ISS_ECI(i,:)' - r_s_eci) / norm(r_ISS_ECI(i,:)' - r_s_eci);
    RA = mod(atan2(u_eci(2), u_eci(1)), 2*pi);
    Dec = asin(u_eci(3));
    RA_Dec_true(i,:) = [RA, Dec];
    RA_Dec_pert(i,:) = [RA, Dec] + sigma * randn(1,2);
end

%% -------------- Printing Tables and Performing IOD/POD ------------------

constants;

% Extract sample epochs and indices at t = 0, 60, 120 seconds (IOD) and every 60 seconds from t = 60 (POD)
idx_IOD = [1, 61, 121];
idx_POD = 121:60:length(t);

% Extract measurement times
date_meas_IOD = date_vec(idx_IOD);
date_meas_POD = date_vec(idx_POD);

% Extract RA/Dec and Az/El (true and perturbed) in radians for IOD
RA_Dec_rad_true_IOD = RA_Dec_true(idx_IOD, :);
RA_Dec_rad_pert_IOD = RA_Dec_pert(idx_IOD, :);
Az_El_rad_true_IOD  = Az_El_true(idx_IOD, :);
Az_El_rad_pert_IOD  = Az_El_pert(idx_IOD, :);

% Extract RA/Dec and Az/El (true and perturbed) in radians for POD
RA_Dec_rad_true_POD = RA_Dec_true(idx_POD, :);
RA_Dec_rad_pert_POD = RA_Dec_pert(idx_POD, :);
Az_El_rad_true_POD  = Az_El_true(idx_POD, :);
Az_El_rad_pert_POD  = Az_El_pert(idx_POD, :);

% Initial orbit estimation using Laplace method (IOD)
x_ISS_ECI_IOD = laplace_orbit_fit(r_s_TH, L, lon, date_meas_IOD, RA_Dec_rad_pert_IOD);

% Initialize for POD estimation
X0_i = x_ISS_ECI_IOD;  % Initial guess
t0   = idx_POD(1);     % Start time for POD

% Gauss-Newton parameters
max_iter  = 100;
epsilon   = 1e-5;
prev_RMS  = 1e3;
iter      = 1;
cond_thresh = 1e11;

% Gauss-Newton Loop
while iter < max_iter
    H = [];  % Jacobian
    r = [];  % Residuals

    for i = 2:length(date_meas_POD)
        tf_POD = idx_POD(i);
        date_POD = date_meas_POD(i);

        % Perturbed Az/El measurement
        Y = Az_El_rad_pert_POD(i, :)';
        Y_hat = compute_AzEl(X0_i, t0, tf_POD, r_s_TH, L, lon, date_POD);
        r_i = Y - Y_hat;

        % Jacobian via finite differencing
        H_i = zeros(2, length(X0_i));
        for k = 1:length(X0_i)
            dx = 0.01 * X0_i(k);
            X_plus  = X0_i; X_plus(k)  = X0_i(k) + dx;
            X_minus = X0_i; X_minus(k) = X0_i(k) - dx;

            az_el_plus  = compute_AzEl(X_plus,  t0, tf_POD, r_s_TH, L, lon, date_POD);
            az_el_minus = compute_AzEl(X_minus, t0, tf_POD, r_s_TH, L, lon, date_POD);

            H_i(:, k) = (az_el_plus - az_el_minus) / (2 * dx);
        end

        H = [H; H_i];
        r = [r; r_i];
    end

    % Gauss-Newton update
    H_t_H = H' * H;         % Compute H' * H
    dX = H_t_H \ (H' * r);  % Gauss-Newton update
    X0_i = X0_i + dX;
    
    cond_num = cond(H_t_H);
    
    if cond_num > cond_thresh
        fprintf('Warning: High condition number detected: %e\n', cond_num);
    end

    RMS = sqrt(r' * r / (2 * length(date_meas_POD)));
    fprintf('Iteration %d, RMS Residual = %.6e\n', iter, RMS);

    RMS_change = abs(RMS - prev_RMS) / prev_RMS;
    if RMS < epsilon || RMS_change < epsilon
        fprintf('Converged at iteration %d.\n', iter);
        break;
    end

    prev_RMS = RMS;
    iter = iter + 1;
end

%% --------------- Convert Radians to Degrees for IOD and POD -------------

vars = {'RA_Dec', 'Az_El'};
types = {'true', 'pert'};
cases = {'IOD', 'POD'};

for v = 1:length(vars)
    for t = 1:length(types)
        for c = 1:length(cases)
            rad_name = sprintf('%s_rad_%s_%s', vars{v}, types{t}, cases{c});
            deg_name = sprintf('%s_deg_%s_%s', vars{v}, types{t}, cases{c});
            eval([deg_name ' = round(rad2deg(' rad_name '), 3);']);
        end
    end
end

%% -------------------------- Create LOS Tables ---------------------------

LOS_IOD_true = table(date_meas_IOD, ...
    RA_Dec_deg_true_IOD(:,1), RA_Dec_deg_true_IOD(:,2), ...
    Az_El_deg_true_IOD(:,1), Az_El_deg_true_IOD(:,2), ...
    'VariableNames', {'Epoch', 'Right Ascension [deg]', 'Declination [deg]', ...
                      'Azimuth [deg]', 'Elevation [deg]'});

LOS_IOD_pert = table(date_meas_IOD, ...
    RA_Dec_deg_pert_IOD(:,1), RA_Dec_deg_pert_IOD(:,2), ...
    Az_El_deg_pert_IOD(:,1), Az_El_deg_pert_IOD(:,2), ...
    'VariableNames', {'Epoch', 'Right Ascension [deg]', 'Declination [deg]', ...
                      'Azimuth [deg]', 'Elevation [deg]'});

LOS_POD_true = table(date_meas_POD, ...
    RA_Dec_deg_true_POD(:,1), RA_Dec_deg_true_POD(:,2), ...
    Az_El_deg_true_POD(:,1), Az_El_deg_true_POD(:,2), ...
    'VariableNames', {'Epoch', 'Right Ascension [deg]', 'Declination [deg]', ...
                      'Azimuth [deg]', 'Elevation [deg]'});

LOS_POD_pert = table(date_meas_POD, ...
    RA_Dec_deg_pert_POD(:,1), RA_Dec_deg_pert_POD(:,2), ...
    Az_El_deg_pert_POD(:,1), Az_El_deg_pert_POD(:,2), ...
    'VariableNames', {'Epoch', 'Right Ascension [deg]', 'Declination [deg]', ...
                      'Azimuth [deg]', 'Elevation [deg]'});

%% ---------------------- State Vector Comparison -------------------------

% True and estimated states
x_ISS_ECI_2   = x_ISS_ECI(61,:)';   % True state
x_ISS_ECI_POD = X0_i;               % POD estimate

% Compute absolute and percentage errors
states = {'IOD', 'POD'};
errors = struct();
for i = 1:length(states)
    est = eval(sprintf('x_ISS_ECI_%s', states{i}));
    abs_err = abs(est - x_ISS_ECI_2);
    pct_err = (abs_err ./ abs(x_ISS_ECI_2)) * 100;
    errors.(sprintf('abs_%s', states{i})) = abs_err;
    errors.(sprintf('pct_%s', states{i})) = pct_err;
end

% Create table
state_components = {'r_x', 'r_y', 'r_z', 'v_x', 'v_y', 'v_z'}';
T = table(state_components, ...
          x_ISS_ECI_IOD, x_ISS_ECI_POD, x_ISS_ECI_2, ...
          errors.pct_IOD, errors.pct_POD, ...
          'VariableNames', {'State Component', 'IOD State', 'POD State', ...
                            'True State', 'IOD State Error (%)', 'POD State Error (%)'});

%% ---------------------- Norm Errors for IOD and POD ---------------------

IOD_position_error_norm = norm(x_ISS_ECI_IOD(1:3) - x_ISS_ECI_2(1:3));
IOD_velocity_error_norm = norm(x_ISS_ECI_IOD(4:6) - x_ISS_ECI_2(4:6));
IOD_total_error_norm    = norm(x_ISS_ECI_IOD - x_ISS_ECI_2);

POD_position_error_norm = norm(x_ISS_ECI_POD(1:3) - x_ISS_ECI_2(1:3));
POD_velocity_error_norm = norm(x_ISS_ECI_POD(4:6) - x_ISS_ECI_2(4:6));
POD_total_error_norm    = norm(x_ISS_ECI_POD - x_ISS_ECI_2);

% Norms of truth data for percent error computation
truth_position_norm = norm(x_ISS_ECI_2(1:3));
truth_velocity_norm = norm(x_ISS_ECI_2(4:6));
truth_total_norm    = norm(x_ISS_ECI_2);

% Percent Errors
IOD_position_percent_error = (IOD_position_error_norm / truth_position_norm) * 100;
IOD_velocity_percent_error = (IOD_velocity_error_norm / truth_velocity_norm) * 100;
IOD_total_percent_error    = (IOD_total_error_norm / truth_total_norm) * 100;

POD_position_percent_error = (POD_position_error_norm / truth_position_norm) * 100;
POD_velocity_percent_error = (POD_velocity_error_norm / truth_velocity_norm) * 100;
POD_total_percent_error    = (POD_total_error_norm / truth_total_norm) * 100;

% Final comparison table with percent error included
Comparison_Error_Table = table(...
    ["Position"; "Velocity"; "Total"], ...
    [IOD_position_error_norm; IOD_velocity_error_norm; IOD_total_error_norm], ...
    [POD_position_error_norm; POD_velocity_error_norm; POD_total_error_norm], ...
    [IOD_position_percent_error; IOD_velocity_percent_error; IOD_total_percent_error], ...
    [POD_position_percent_error; POD_velocity_percent_error; POD_total_percent_error], ...
    {'km'; 'km/s'; 'mixed'}, ...
    'VariableNames', {'Component', 'IOD Error', 'POD Error', ...
                      'IOD % Error', 'POD % Error', 'Units'});

%% ------------------------ Display Results -------------------------------

disp('--- Random Seed Number ---');               disp(rng().Seed);
disp('--- TRUE LOS MEASUREMENTS (IOD) ---');      disp(LOS_IOD_true);
disp('--- PERTURBED LOS MEASUREMENTS (IOD) ---'); disp(LOS_IOD_pert);
disp('--- TRUE LOS MEASUREMENTS (POD) ---');      disp(LOS_POD_true);
disp('--- PERTURBED LOS MEASUREMENTS (POD) ---'); disp(LOS_POD_pert);

disp('--- COMPARISON OF IOD AND POD STATES ---');
disp(T);

disp('--- COMPARISON OF IOD AND POD STATE ESTIMATION ERRORS ---');
disp(Comparison_Error_Table);

%% Plotting

t_prop = [0 24*60*60];

[t, x_ISS_ECI_prop] = ode45(@(t,y) two_body_ode(t, y), t_prop, x_ISS_ECI_2, opts);
[~, x_ISS_ECI_IOD_prop] = ode45(@(t,y) two_body_ode(t, y), t, x_ISS_ECI_IOD, opts);
[~, x_ISS_ECI_POD_prop] = ode45(@(t,y) two_body_ode(t, y), t, x_ISS_ECI_POD, opts);

IOD_error = x_ISS_ECI_prop - x_ISS_ECI_IOD_prop;
POD_error = x_ISS_ECI_prop - x_ISS_ECI_POD_prop;

% Plotting state errors in a 6x2 subplot

fs = 16;
labels = {'x [km]', 'y [km]', 'z [km]', 'vx [km/s]', 'vy [km/s]', 'vz [km/s]'};

figure;
for i = 1:6
    % Left column: IOD error
    subplot(6,2,2*i-1);
    plot(t/3600, IOD_error(:,i), 'k', 'LineWidth', 2); % black line
    grid on;
    if i == 1
        title('IOD', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
    end
    if i == 6
        xlabel('Time [hr]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
    end
    ylabel(labels{i}, 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
    set(gca, 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');

    % Right column: POD error
    subplot(6,2,2*i);
    plot(t/3600, POD_error(:,i), 'k', 'LineWidth', 2); % black line
    grid on;
    if i == 1
        title('POD', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
    end
    if i == 6
        xlabel('Time [hr]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
    end
    % No ylabel for right column
    set(gca, 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
end

sgtitle('IOD and POD State Error Comparison', ...
    'FontName', 'Times New Roman', 'FontSize', fs + 2, 'FontWeight', 'bold');