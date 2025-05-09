clear
close all
clc

% Populating constants
constants

% Set font size
fs = 14;  % You can change this to any size you prefer

% Set line width
lineWidth = 1.5;

% Load the data from the CSV file
data = readmatrix('ME594 centroid data for Final.xlsx');

% Extracting coordinates for each image
x_cam_1 = data(:, 1); % Image 1 x_cam
y_cam_1 = data(:, 2); % Image 1 y_cam

x_cam_2 = data(:, 4); % Image 2 x_cam
y_cam_2 = data(:, 5); % Image 2 y_cam

x_cam_3 = data(:, 7); % Image 3 x_cam
y_cam_3 = data(:, 8); % Image 3 y_cam

x_cam_4 = data(:, 10); % Image 4 x_cam
y_cam_4 = data(:, 11); % Image 4 y_cam

x_cam_5 = data(:, 13); % Image 5 x_cam
y_cam_5 = data(:, 14); % Image 5 y_cam

% Plotting the coordinates for each camera set with different colors
figure;
hold on;

% Plot Image 1 in red with linewidth 1.5
plot(x_cam_1, y_cam_1, 'ro', 'DisplayName', 'Image 1', 'LineWidth', lineWidth);

% Plot Image 2 in blue with linewidth 1.5
plot(x_cam_2, y_cam_2, 'bo', 'DisplayName', 'Image 2', 'LineWidth', lineWidth);

% Plot Image 3 in green with linewidth 1.5
plot(x_cam_3, y_cam_3, 'go', 'DisplayName', 'Image 3', 'LineWidth', lineWidth);

% Plot Image 4 in purple with linewidth 1.5
plot(x_cam_4, y_cam_4, 'mo', 'DisplayName', 'Image 4', 'LineWidth', lineWidth);

% Plot Image 5 in orange with linewidth 1.5
plot(x_cam_5, y_cam_5, 'co', 'DisplayName', 'Image 5', 'LineWidth', lineWidth);

% Add title and labels
title('Camera Coordinates for Images 1-5', 'FontSize', fs, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
xlabel('X Coordinates (pixels)', 'FontSize', fs, 'FontWeight', 'bold', 'FontName', 'Times New Roman');
ylabel('Y Coordinates (pixels)', 'FontSize', fs, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

% Add legend
legend show;

% Enable grid
grid on;

% Set axis properties
set(gca, 'FontSize', fs, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

% Hold off to stop adding to the current plot
hold off;

%% Problem 2

% rng('shuffle');
rng(464876754);
currentState = rng;
disp(currentState.Seed);

% Starting and ending dates
date_0 = datetime('2025-05-05 03:18:34', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');

% Latitude longitude for Albuquerque, New Mexico
L = deg2rad(35.0844);   
lon = deg2rad(-106.6504);

% ISS Initial State
r_0_ISS_ECI = [5165.8127; 
              -1938.2678;
              -3951.6672];  % Position [km]
v_0_ISS_ECI = [4.7525;
               4.4652;
               4.0249];     % Velocity [km/s]
x_0_ISS_ECI = [r_0_ISS_ECI; v_0_ISS_ECI]; % [km; km/s]

% Specifying time since epoch to numerically integrate
t_span = 0:60*2;%

% Integrate the ODEs
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t, x_ISS_ECI] = ode45(@(t, y) two_body_ode(t, y), t_span, x_0_ISS_ECI, options);

r_ISS_ECI = x_ISS_ECI(:,1:3); %Extracting ISS position

N = length(t);                          % Number of measurements
RA_Dec_true = zeros(N,2);               % Allocate matrix for RA/Dec (true)
Az_El_true  = zeros(N,2);               % Allocate matrix for Az/El (true)
RA_Dec_pert = zeros(N,2);               % Allocate matrix for RA/Dec (perturbed)
Az_El_pert  = zeros(N,2);               % Allocate matrix for Az/El (perturbed)
sigma = deg2rad(0.01);                  % Measurement noise standard deviation [rad]
date_vec = date_0 + seconds(t);         % Create datetime array

for i = 1:N
    date = date_0 + seconds(t(i));      % Current date
    theta_g = GMST(date);               % GMST
    theta = lon + theta_g;

    Q_TH_ECI = rot_ijk_sez(L, theta);   % ECI->TH rotation matrix

    r_ISS_TH = Q_TH_ECI * r_ISS_ECI(i,:)';      % ISS position in TH frame
    r_s_ISS_TH = r_ISS_TH - r_s_TH;             % Vector from sensor to ISS
    u_s_ISS_TH = r_s_ISS_TH / norm(r_s_ISS_TH); % Unit vector to ISS in TH

    % Azimuth and Elevation from TH frame
    el = asin(-u_s_ISS_TH(3));                % Elevation: arcsin(-z)
    az = atan2(u_s_ISS_TH(1), u_s_ISS_TH(2)); % Azimuth: atan2(x, y)
    if az < 0
        az = az + 2*pi;
    end
    Az_El_true(i,:) = [az, el];               % Store true Az/El

    % Perturbed Az/El
    Az_El_pert(i,1) = az + sigma * randn;     % Perturbed Azimuth
    Az_El_pert(i,2) = el + sigma * randn;     % Perturbed Elevation

    % Sensor position in ECI
    r_s_ECI = Q_TH_ECI' * r_s_TH;
    r_s_ISS_ECI = r_ISS_ECI(i,:)' - r_s_ECI;
    u_s_ISS_ECI = r_s_ISS_ECI / norm(r_s_ISS_ECI); % Unit vector to ISS in ECI frame

    % RA/Dec from ECI
    RA  = atan2(u_s_ISS_ECI(2), u_s_ISS_ECI(1));
    if RA < 0
        RA = RA + 2*pi;
    end
    Dec = asin(u_s_ISS_ECI(3));

    RA_Dec_true(i,:) = [RA, Dec];             % Store true RA/Dec

    % Perturbed RA/Dec
    RA_Dec_pert(i,1) = RA  + sigma * randn;   % Perturbed RA
    RA_Dec_pert(i,2) = Dec + sigma * randn;   % Perturbed Dec
end

%% Printing Tables

constants;
% Extract sample epochs and indices at t = 0, 60, 120 seconds (IOD) and every 3 seconds from t=60 (POD)
idx_IOD = [1, 61, length(t)];
idx_POD = 61:3:length(t);

% Extract epoch time vectors
date_meas_IOD = date_vec(idx_IOD);
date_meas_POD = date_vec(idx_POD);
RA_Dec_pert_IOD = RA_Dec_pert(idx_IOD, :);

% Extract RA/Dec and Az/El (true and perturbed) in radians for IOD
RA_Dec_rad_true_IOD  = RA_Dec_true(idx_IOD, :);
Az_El_rad_true_IOD   = Az_El_true(idx_IOD, :);
RA_Dec_rad_pert_IOD  = RA_Dec_pert(idx_IOD, :);
Az_El_rad_pert_IOD   = Az_El_pert(idx_IOD, :);

% Extract RA/Dec and Az/El (true and perturbed) in radians for POD
RA_Dec_rad_true_POD  = RA_Dec_true(idx_POD, :);
Az_El_rad_true_POD   = Az_El_true(idx_POD, :);
RA_Dec_rad_pert_POD  = RA_Dec_pert(idx_POD, :);
Az_El_rad_pert_POD   = Az_El_pert(idx_POD, :);

x_ISS_ECI_IOD = laplace_orbit_fit(r_s_TH, L, lon, date_meas_IOD, RA_Dec_pert_IOD);
X0 = x_ISS_ECI_IOD;
X0_i = X0;
t0 = idx_POD(1);

epsilon = 1e-6;

for i = 2:length(date_meas_POD)
    tf_POD = idx_POD(i);
    date_POD = date_meas_POD(i);

    Y     = Az_El_rad_pert_POD(i,:)';
    Y_hat = compute_AzEl(X0_i, t0, tf_POD, r_s_TH, L, lon, date_POD);
    r     = Y - Y_hat;

    H_0 = zeros(2, 6);
    
    for k = 1:length(X0_i)
        dx_pert_i = X0_i(k) * 0.01;

        X0_pert_plus = X0_i;
        X0_pert_plus(k) = X0_i(k) + dx_pert_i;

        X0_pert_minus = X0_i;
        X0_pert_minus(k) = X0_i(k) - dx_pert_i;

        az_el_plus = compute_AzEl(X0_pert_plus, t0, tf_POD, r_s_TH, L, lon, date_POD);
        az_el_minus = compute_AzEl(X0_pert_minus, t0, tf_POD, r_s_TH, L, lon, date_POD);

        H_0(:, k) = (az_el_plus - az_el_minus) / (2 * dx_pert_i);
    end

    cond_num = cond(H_0' * H_0);
    criterion = norm(r);

    dX_i = pinv(H_0' * H_0) * H_0' * r;
    X0_i = X0_i + dX_i;

    if  criterion < epsilon
        fprintf('Breaking at iteration %d because residual norm is below tolerance.\n', i);
        break;
    end
end

% Convert extracted radian values to degrees and round
RA_Dec_deg_true_IOD = round(rad2deg(RA_Dec_rad_true_IOD), 3);
Az_El_deg_true_IOD  = round(rad2deg(Az_El_rad_true_IOD), 3);
RA_Dec_deg_pert_IOD = round(rad2deg(RA_Dec_rad_pert_IOD), 3);
Az_El_deg_pert_IOD  = round(rad2deg(Az_El_rad_pert_IOD), 3);

RA_Dec_deg_true_POD = round(rad2deg(RA_Dec_rad_true_POD), 3);
Az_El_deg_true_POD  = round(rad2deg(Az_El_rad_true_POD), 3);
RA_Dec_deg_pert_POD = round(rad2deg(RA_Dec_rad_pert_POD), 3);
Az_El_deg_pert_POD  = round(rad2deg(Az_El_rad_pert_POD), 3);

% Create IOD tables
LOS_IOD_true = table(date_meas_IOD, ...
    RA_Dec_deg_true_IOD(:,1), RA_Dec_deg_true_IOD(:,2), ...
    Az_El_deg_true_IOD(:,1), Az_El_deg_true_IOD(:,2), ...
    'VariableNames', {'Epoch', ...
                      'Right Ascension [deg]', 'Declination [deg]', ...
                      'Azimuth [deg]', 'Elevation [deg]'});

LOS_IOD_pert = table(date_meas_IOD, ...
    RA_Dec_deg_pert_IOD(:,1), RA_Dec_deg_pert_IOD(:,2), ...
    Az_El_deg_pert_IOD(:,1), Az_El_deg_pert_IOD(:,2), ...
    'VariableNames', {'Epoch', ...
                      'Right Ascension [deg]', 'Declination [deg]', ...
                      'Azimuth [deg]', 'Elevation [deg]'});

% Create POD tables
LOS_POD_true = table(date_meas_POD, ...
    RA_Dec_deg_true_POD(:,1), RA_Dec_deg_true_POD(:,2), ...
    Az_El_deg_true_POD(:,1), Az_El_deg_true_POD(:,2), ...
    'VariableNames', {'Epoch', ...
                      'Right Ascension [deg]', 'Declination [deg]', ...
                      'Azimuth [deg]', 'Elevation [deg]'});

LOS_POD_pert = table(date_meas_POD, ...
    RA_Dec_deg_pert_POD(:,1), RA_Dec_deg_pert_POD(:,2), ...
    Az_El_deg_pert_POD(:,1), Az_El_deg_pert_POD(:,2), ...
    'VariableNames', {'Epoch', ...
                      'Right Ascension [deg]', 'Declination [deg]', ...
                      'Azimuth [deg]', 'Elevation [deg]'});

x_ISS_ECI_2 = x_ISS_ECI(61,:)';

% Assign POD state from the previous estimation
x_ISS_ECI_POD = X0_i;  % POD state values after estimation

% Calculate the error (absolute error) for IOD, POD, and true state
error_IOD = abs(x_ISS_ECI_IOD - x_ISS_ECI_2);
error_POD = abs(x_ISS_ECI_POD - x_ISS_ECI_2);

% Calculate the percent error for IOD and POD states
percent_error_IOD = (error_IOD ./ abs(x_ISS_ECI_2)) * 100;
percent_error_POD = (error_POD ./ abs(x_ISS_ECI_2)) * 100;

% Create a table with the values for IOD, POD, and True state comparison
state_components = {'r_x', 'r_y', 'r_z', 'v_x', 'v_y', 'v_z'};
IOD_state = x_ISS_ECI_IOD;     % OD state values (IOD)
True_state = x_ISS_ECI_2;     % True state values
POD_state = x_ISS_ECI_POD;    % POD state values

State_error_IOD = percent_error_IOD;  % Percent error for IOD
State_error_POD = percent_error_POD; % Percent error for POD

% Display the table with the state components and errors
T = table(state_components', IOD_state, True_state, State_error_IOD, State_error_POD, ...
    'VariableNames', {'State Component', 'IOD State', 'True State', 'IOD State Error (%)', 'POD State Error (%)'});

% Display tables
disp('--- TRUE LOS MEASUREMENTS (IOD) ---');
disp(LOS_IOD_true);

disp('--- PERTURBED LOS MEASUREMENTS (IOD) ---');
disp(LOS_IOD_pert);

disp('--- TRUE LOS MEASUREMENTS (POD) ---');
disp(LOS_POD_true);

disp('--- PERTURBED LOS MEASUREMENTS (POD) ---');
disp(LOS_POD_pert);

% Display the table
disp(T);