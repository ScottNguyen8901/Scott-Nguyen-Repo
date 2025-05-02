clear
close all
clc

% Populating constants
constants

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
    u_s_ISS_ECI = r_s_ISS_ECI / norm(r_s_ISS_ECI);

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

% Extract sample epochs and indices at t = 0, 60, 120 seconds
idx_sample = [1, 61, length(t)];
date_meas = date_vec(idx_sample);

% Extract RA/Dec and Az/El (true and perturbed)
RA_Dec_meas_true  = RA_Dec_true(idx_sample, :);
Az_El_meas_true   = Az_El_true(idx_sample, :);
RA_Dec_meas_pert  = RA_Dec_pert(idx_sample, :);
Az_El_meas_pert   = Az_El_pert(idx_sample, :);

% Convert from radians to degrees and round
RA_Dec_deg_true  = round(rad2deg(RA_Dec_meas_true), 3);
Az_El_deg_true   = round(rad2deg(Az_El_meas_true), 3);
RA_Dec_deg_pert  = round(rad2deg(RA_Dec_meas_pert), 3);
Az_El_deg_pert   = round(rad2deg(Az_El_meas_pert), 3);

% Create tables
LOS_table_true = table(date_meas, ...
    RA_Dec_deg_true(:,1), RA_Dec_deg_true(:,2), ...
    Az_El_deg_true(:,1), Az_El_deg_true(:,2), ...
    'VariableNames', {'Epoch', ...
                      'Right Ascension [deg]', 'Declination [deg]', ...
                      'Azimuth [deg]', 'Elevation [deg]'});

LOS_table_pert = table(date_meas, ...
    RA_Dec_deg_pert(:,1), RA_Dec_deg_pert(:,2), ...
    Az_El_deg_pert(:,1), Az_El_deg_pert(:,2), ...
    'VariableNames', {'Epoch', ...
                      'Right Ascension [deg]', 'Declination [deg]', ...
                      'Azimuth [deg]', 'Elevation [deg]'});

% Display both tables
disp('--- TRUE LOS MEASUREMENTS ---');
disp(LOS_table_true);

disp('--- PERTURBED LOS MEASUREMENTS ---');
disp(LOS_table_pert);

x_ISS_ECI_OD = laplace_orbit_fit(r_s_TH, L, lon, date_meas, RA_Dec_meas_pert);
x_ISS_ECI_2 = x_ISS_ECI(61,:)';

% Calculate the error (absolute error)
error = abs(x_ISS_ECI_OD - x_ISS_ECI_2);

% Calculate the percent error
percent_error = (error ./ abs(x_ISS_ECI_2)) * 100;

% Create a table with the values
state_components = {'r_x', 'r_y', 'r_z', 'v_x', 'v_y', 'v_z'};
OD_state = x_ISS_ECI_OD;      % OD state values
True_state = x_ISS_ECI_2;     % True state values
State_error = percent_error;  % Percent error values

% Display the table
T = table(state_components', OD_state, True_state, State_error, ...
    'VariableNames', {'State Component', 'OD State', 'True State', 'State Error (%)'});

% Display the table
disp(T);

%% Propogating further in time

P = 2*pi * sqrt((6785.6^3) / mu_E);  % Orbital period

t_span_prop = [0 5*P];
[t_prop, x_ISS_ECI_OD_Prop] = ...
    ode45(@(t, y) two_body_ode(t, y), t_span_prop, x_ISS_ECI_OD, options);
[~, x_ISS_ECI_OD_True] = ...
    ode45(@(t, y) two_body_ode(t, y), t_prop, x_ISS_ECI_2, options);

error_prop = x_ISS_ECI_OD_True - x_ISS_ECI_OD_Prop;

%% Plotting

% Define font size variable
fs = 24;

% Create the 3x2 subplot
figure;

% Plot the error in the x-component (position)
subplot(3, 2, 1);
plot(t_prop / P, error_prop(:, 1), 'k', 'LineWidth', 2);
title('e_{r_x} [km]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
xlabel('Time [P]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
ylabel('Error [km]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
grid on;

% Plot the error in the y-component (position)
subplot(3, 2, 3);
plot(t_prop / P, error_prop(:, 2), 'k', 'LineWidth', 2);
title('e_{r_y} [km]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
xlabel('Time [P]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
ylabel('Error [km]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
grid on;

% Plot the error in the z-component (position)
subplot(3, 2, 5);
plot(t_prop / P, error_prop(:, 3), 'k', 'LineWidth', 2);
title('e_{r_z} [km]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
xlabel('Time [P]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
ylabel('Error [km]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
grid on;

% Plot the error in the vx-component (velocity)
subplot(3, 2, 2);
plot(t_prop / P, error_prop(:, 4), 'k', 'LineWidth', 2);
title('e_{v_x} [km/s]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
xlabel('Time [P]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
ylabel('Error [km/s]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
grid on;

% Plot the error in the vy-component (velocity)
subplot(3, 2, 4);
plot(t_prop / P, error_prop(:, 5), 'k', 'LineWidth', 2);
title('e_{v_y} [km/s]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
xlabel('Time [P]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
ylabel('Error [km/s]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
grid on;

% Plot the error in the vz-component (velocity)
subplot(3, 2, 6);
plot(t_prop / P, error_prop(:, 6), 'k', 'LineWidth', 2);
title('e_{v_z} [km/s]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
xlabel('Time [P]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
ylabel('Error [km/s]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
grid on;

sgtitle('Perturbed State Error', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', 32);