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
r_s_TH = [0; 0; R_E]; % Location of sensor in TH (SEZ) frame      

r_0_ISS_ECI = [5165.8127; 
              -1938.2678;
              -3951.6672];  % Position in km
v_0_ISS_ECI = [4.7525;
               4.4652;
               4.0249];  % Velocity in km/s
x_0_ISS_ECI = [r_0_ISS_ECI; v_0_ISS_ECI];

% Specifying time since epoch to numerically integrate
t_span = 0:60*2;%

% Integrate the ODEs
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t, x_ISS_ECI] = ode45(@(t, y) two_body_ode(t, y, mu_E), t_span, x_0_ISS_ECI, options);

r_ISS_ECI = x_ISS_ECI(:,1:3); %Extracting ISS position

N = length(t);         % Number of measurements
RA_Dec = zeros(N,2); % Allocating matrix for LOS measurements
Az_El = zeros(N,2);    % Allocating matrix for Az/El measurements
date_vec = date_0 + seconds(t); % Create an array of datetime objects based on time span

for i = 1:N
    date = date_0 + seconds(t(i));                     % Current date
    theta_g = GMST(date);                              % GMST
    theta = lon + theta_g;

    Q_TH_ECI = rot_ijk_sez(L, theta);                  % ECI->TH

    r_ISS_TH = Q_TH_ECI * r_ISS_ECI(i,:)';             % ISS position TH
    r_s_ISS_TH = r_ISS_TH - r_s_TH;                    % Sens to ISS in TH
    u_s_ISS_TH = r_s_ISS_TH / norm(r_s_ISS_TH);        % Sens to ISS unit vec in TH

    % Azimuth and Elevation from TH frame
    el = asin(-u_s_ISS_TH(3));                         % Elevation: arcsin(-z)
    az = atan2(u_s_ISS_TH(1), u_s_ISS_TH(2));          % Azimuth: atan2(x, y)
    if az < 0
        az = az + 2*pi;                                % Ensure Azimuth is in [0, 2pi]
    end
    Az_El(i,:) = [az, el];                         % Store Azimuth and Elevation

    % Convert r_s_TH back to ECI to compute RA/Dec
    r_s_ECI = Q_TH_ECI' * r_s_TH;                      % Position of sensor in ECI
    r_s_ISS_ECI = r_ISS_ECI(i,:)' - r_s_ECI;           % Vector from sensor to ISS in ECI
    u_s_ISS_ECI = r_s_ISS_ECI / norm(r_s_ISS_ECI);     % Unit vector in ECI

    % Compute RA and Dec from ECI vector
    RA  = atan2(u_s_ISS_ECI(2), u_s_ISS_ECI(1));        % RA: arctangent of y/x
    if RA < 0
        RA = RA + 2*pi;                                 % Ensure RA is in [0, 2pi]
    end
    Dec = asin(u_s_ISS_ECI(3));                         % Declination: arcsin(z)

    RA_Dec(i,:) = [RA, Dec];                            % Store in output matrix
end

%% Plotting

% Extracting epoch and LOS measurements at t = 0, 60, 120 s
date_meas = [date_vec(1);
             date_vec(61);
             date_vec(end)];

RA_Dec_meas = [RA_Dec(1,:);
               RA_Dec(61,:);
               RA_Dec(end,:)];

Az_El_meas = [Az_El(1,:);
              Az_El(61,:);
              Az_El(end,:)];

% Convert LOS measurements from radians to degrees and round to 3 decimals
RA_Dec_deg = round(rad2deg(RA_Dec_meas), 3);
Az_El_deg = round(rad2deg(Az_El_meas), 3);

% Create a combined table with appropriate column names
LOS_table = table(date_meas, ...
    RA_Dec_deg(:,1), RA_Dec_deg(:,2), ...
    Az_El_deg(:,1), Az_El_deg(:,2), ...
    'VariableNames', {'Epoch', ...
                      'Right Ascension [deg]', 'Declination [deg]', ...
                      'Azimuth [deg]', 'Elevation [deg]'});

% Display the table
disp(LOS_table);

x_ISS_ECI_OD = laplace_orbit_fit(r_s_TH, L, lon, date_meas, RA_Dec_meas);
x_ISS_ECI_2 = x_ISS_ECI(61,:)';

% Calculate the error (absolute error)
error = abs(x_ISS_ECI_OD - x_ISS_ECI_2);

% Calculate the percent error
percent_error = (error ./ abs(x_ISS_ECI_2)) * 100;

% Create a table with the values
state_components = {'r_x', 'r_y', 'r_z', 'v_x', 'v_y', 'v_z'};  % State component names
OD_state = x_ISS_ECI_OD;  % OD state values
True_state = x_ISS_ECI_2;  % True state values
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
    ode45(@(t, y) two_body_ode(t, y, mu_E), t_span_prop, x_ISS_ECI_OD, options);
[~, x_ISS_ECI_OD_True] = ...
    ode45(@(t, y) two_body_ode(t, y, mu_E), t_prop, x_ISS_ECI_2, options);

error_prop = x_ISS_ECI_OD_True - x_ISS_ECI_OD_Prop;

%% Plotting

% Define font size variable
fs = 24;

% Create the 3x2 subplot
figure;

% Plot the error in the x-component (position)
subplot(3, 2, 1);
plot(t_prop / P, error_prop(:, 1), 'k', 'LineWidth', 2);  % Solid black line
title('e_{r_x} [km]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
xlabel('Time [P]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
ylabel('Error [km]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
grid on;

% Plot the error in the y-component (position)
subplot(3, 2, 3);
plot(t_prop / P, error_prop(:, 2), 'k', 'LineWidth', 2);  % Solid black line
title('e_{r_y} [km]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
xlabel('Time [P]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
ylabel('Error [km]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
grid on;

% Plot the error in the z-component (position)
subplot(3, 2, 5);
plot(t_prop / P, error_prop(:, 3), 'k', 'LineWidth', 2);  % Solid black line
title('e_{r_z} [km]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
xlabel('Time [P]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
ylabel('Error [km]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
grid on;

% Plot the error in the vx-component (velocity)
subplot(3, 2, 2);
plot(t_prop / P, error_prop(:, 4), 'k', 'LineWidth', 2);  % Solid black line
title('e_{v_x} [km/s]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
xlabel('Time [P]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
ylabel('Error [km/s]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
grid on;

% Plot the error in the vy-component (velocity)
subplot(3, 2, 4);
plot(t_prop / P, error_prop(:, 5), 'k', 'LineWidth', 2);  % Solid black line
title('e_{v_y} [km/s]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
xlabel('Time [P]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
ylabel('Error [km/s]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
grid on;

% Plot the error in the vz-component (velocity)
subplot(3, 2, 6);
plot(t_prop / P, error_prop(:, 6), 'k', 'LineWidth', 2);  % Solid black line
title('e_{v_z} [km/s]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
xlabel('Time [P]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
ylabel('Error [km/s]', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', fs);
grid on;

sgtitle('State Error', 'FontName', 'Times New Roman', 'FontWeight', 'bold', 'FontSize', 32);