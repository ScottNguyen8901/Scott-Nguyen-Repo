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
t_span = 0:60*2;

% Integrate the ODEs
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t, x_ISS_ECI] =...
    ode45(@(t, y) two_body_ode(t, y), t_span, x_0_ISS_ECI, options);

r_ISS_ECI = x_ISS_ECI(:,1:3); % Extracting ISS position

num_runs = 1000;                        % Number of Monte Carlo simulations
N = length(t);                          % Number of measurements
RA_Dec_true = zeros(N,2);               % Allocate matrix for RA/Dec (true)
Az_El_true  = zeros(N,2);               % Allocate matrix for Az/El (true)
RA_Dec_pert = zeros(N,2);               % Allocate matrix for RA/Dec (perturbed)
Az_El_pert  = zeros(N,2);               % Allocate matrix for Az/El (perturbed)
sigma = deg2rad(0.01);                  % Measurement noise standard deviation [rad]
date_vec = date_0 + seconds(t);         % Create datetime array

% Initialize array to store error for each Monte Carlo run
state_errors = zeros(num_runs, 6);

for run = 1:num_runs
    % Initialize true and perturbed values for this run
    for j = 1:N
        date = date_0 + seconds(t(j));      % Current date
        theta_g = GMST(date);               % GMST
        theta = lon + theta_g;

        Q_TH_ECI = rot_ijk_sez(L, theta);   % ECI->TH rotation matrix

        r_ISS_TH = Q_TH_ECI * r_ISS_ECI(j,:)';     % ISS position in TH frame
        r_s_ISS_TH = r_ISS_TH - r_s_TH;            % Vector from sensor to ISS
        u_s_ISS_TH = r_s_ISS_TH / norm(r_s_ISS_TH);% Unit vector to ISS in TH

        % Azimuth and Elevation from TH frame
        el = asin(-u_s_ISS_TH(3));                % Elevation: arcsin(-z)
        az = atan2(u_s_ISS_TH(1), u_s_ISS_TH(2)); % Azimuth: atan2(x, y)
        if az < 0
            az = az + 2*pi;
        end
        Az_El_true(j,:) = [az, el];               % Store true Az/El

        % Perturbed Az/El
        Az_El_pert(j,1) = az + sigma * randn;     % Perturbed Azimuth
        Az_El_pert(j,2) = el + sigma * randn;     % Perturbed Elevation

        % Sensor position in ECI
        r_s_ECI = Q_TH_ECI' * r_s_TH;
        r_s_ISS_ECI = r_ISS_ECI(j,:)' - r_s_ECI;
        u_s_ISS_ECI = r_s_ISS_ECI / norm(r_s_ISS_ECI);

        % RA/Dec from ECI
        RA  = atan2(u_s_ISS_ECI(2), u_s_ISS_ECI(1));
        if RA < 0
            RA = RA + 2*pi;
        end
        Dec = asin(u_s_ISS_ECI(3));

        RA_Dec_true(j,:) = [RA, Dec];             % Store true RA/Dec

        % Perturbed RA/Dec
        RA_Dec_pert(j,1) = RA  + sigma * randn;   % Perturbed RA
        RA_Dec_pert(j,2) = Dec + sigma * randn;   % Perturbed Dec
    end

    % Extract sample epochs and indices at t = 0, 60, 120 seconds
    idx_sample = [1, 61, length(t)];
    date_meas = date_vec(idx_sample);
    RA_Dec_meas_pert  = RA_Dec_pert(idx_sample, :);

    % Perform Laplace orbit fit to get the state
    x_ISS_ECI_OD = laplace_orbit_fit(r_s_TH, L, lon, date_meas, RA_Dec_meas_pert);
    x_ISS_ECI_2 = x_ISS_ECI(61,:)' ;

    % Calculate the error (absolute error)
    error = abs(x_ISS_ECI_OD - x_ISS_ECI_2);

    % Store the error for this run
    state_errors(run, :) = error;

    % Print progress
    clc
    fprintf('Progress: %.2f%%\n', (run / num_runs) * 100);
end

%% Plotting

% Define font size variable
fs = 24;

% Plot histogram of state errors
figure;

subplot(3,2,1);
histogram(state_errors(:,1), 20);
title('Position Error (x)', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
xlabel('r_x [km]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
ylabel('Frequency', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
grid on;

subplot(3,2,3);
histogram(state_errors(:,2), 20);
title('Position Error (y)', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
xlabel('r_y [km]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
ylabel('Frequency', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
grid on;

subplot(3,2,5);
histogram(state_errors(:,3), 20);
title('Position Error (z)', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
xlabel('r_z [km]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
ylabel('Frequency', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
grid on;

subplot(3,2,2);
histogram(state_errors(:,4), 20);
title('Velocity Error (v_x)', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
xlabel('v_x [km/s]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
ylabel('Frequency', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
grid on;

subplot(3,2,4);
histogram(state_errors(:,5), 20);
title('Velocity Error (v_y)', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
xlabel('v_y [km/s]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
ylabel('Frequency', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
grid on;

subplot(3,2,6);
histogram(state_errors(:,6), 20);
title('Velocity Error (v_z)', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
xlabel('v_z [km/s]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
ylabel('Frequency', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
grid on;

sgtitle('State Error Histogram', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', 34);