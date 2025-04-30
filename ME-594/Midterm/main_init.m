clear
close all
clc

% Populating constants
constants

% Define the orbital elements in km, radians, seconds
koe.a = 6785.6;              % Semi-major axis [km]
koe.e = 0.00023;             % Eccentricity []
koe.i = deg2rad(51.6367);    % Inclination [rad]
koe.W = deg2rad(13.9687);    % RAAN [rad]
koe.w = deg2rad(82.7408);    % Argument of Periapsis [rad]
koe.f = deg2rad(281.9577);   % True anomaly [rad]

% Starting and ending dates
date_0 = datetime('2025-04-25 17:56:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
date_f = datetime('2025-05-06 00:00:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');

% Latitude longitude for Albuqurque, New Mexico
L = deg2rad(35.0844);   
lon = deg2rad(-106.6504);
r_s_TH = [0; 0; R_E]; % Location of sensor in TH (SEZ) frame      

% Compute x, y, z on the sphere
x = R_E * cos(L) * cos(lon);
y = R_E * cos(L) * sin(lon);
z = R_E * sin(L);

% Converting orbital elements to state vector
[r_0_ISS_ECI, v_0_ISS_ECI] = koe_to_rv(koe, mu_E);
x_0_ISS_ECI = [r_0_ISS_ECI; v_0_ISS_ECI];

% Specifying time since epoch to numerically integrate
tf = seconds(date_f - date_0);
t_span = linspace(0, tf, 1000);%

% Integrate the ODEs
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t, x_ISS_ECI] = ode45(@(t, y) two_body_ode(t, y, mu_E), t_span, x_0_ISS_ECI, options);

r_ISS_ECI = x_ISS_ECI(:,1:3); %Extracting ISS position

N = length(t);             % Number of measurements
LOS_meas = zeros(N,2);     % Allocating matrix for RA/Dec measurements
AzEl_meas = zeros(N,2);    % Allocating matrix for Az/El measurements

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
    AzEl_meas(i,:) = [az, el];                         % Store Azimuth and Elevation

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

    LOS_meas(i,:) = [RA, Dec];                         % Store RA/Dec
end

% Find the maximum elevation and its corresponding time
[max_el, max_el_idx] = max(AzEl_meas(:,2));  % max elevation and its index
t_max_el = t(max_el_idx);                % time at which max elevation occurs

% Display result
fprintf('Maximum Elevation: %.6f radians (%.3f deg)\n', max_el, rad2deg(max_el));
fprintf('Time since Epoch at Max Elevation: %.3f seconds\n', t_max_el);

% Extract the corresponding state from x_ISS_ECI
state_at_dt = x_ISS_ECI(max_el_idx, :);  % State vector at time dt

% Extract position and velocity components
position = state_at_dt(1:3);  % Position (x, y, z) in km
velocity = state_at_dt(4:6);  % Velocity (vx, vy, vz) in km/s

% Extract azimuth and elevation at dt
RA_Dec_at_dt = LOS_meas(max_el_idx, :);  % [Azimuth, Elevation] at dt
RA = RA_Dec_at_dt(1);             % Right Ascension in degrees
Dec = RA_Dec_at_dt(2);            % Declination in degrees

% Calculate the corresponding date for dt
date_at_dt = date_0 + seconds(t_max_el);  % Add seconds to the starting date

% Display the state vector in a neat format
disp('State at dt:');
disp(' ');  % Adding a space for better readability

% Display position
fprintf('Position (x, y, z) [km]:\n');
fprintf('x = %.4f km\n', position(1));
fprintf('y = %.4f km\n', position(2));
fprintf('z = %.4f km\n', position(3));

disp(' ');  % Adding a space for better readability

% Display velocity
fprintf('Velocity (vx, vy, vz) [km/s]:\n');
fprintf('vx = %.4f km/s\n', velocity(1));
fprintf('vy = %.4f km/s\n', velocity(2));
fprintf('vz = %.4f km/s\n', velocity(3));

disp(' ');  % Adding a space for better readability

% Display azimuth and elevation
fprintf('Right Ascension and Declination [deg]:\n');
fprintf('Right Ascension = %.4f rad\n', RA);
fprintf('Declination = %.4f rad\n', Dec);

disp(' ');  % Adding a space for better readability

% Display the corresponding date
fprintf('Corresponding date: %s\n', datestr(date_at_dt));

%% Plotting

% Set the desired font size
fs = 24;  % Change this value to your preferred font size

% Create time vector as datetime
time_vec = date_0 + seconds(t); % datetimes matching each measurement

% Convert Right Ascension and Declination to degrees
RA_deg = rad2deg(LOS_meas(:,1));
Dec_deg = rad2deg(LOS_meas(:,2));

% Convert Azimuth and Elevation to degrees
Az_deg = rad2deg(AzEl_meas(:,1));
El_deg = rad2deg(AzEl_meas(:,2));

% Create grids for interpolation
time_grid = linspace(0, seconds(time_vec(end) - time_vec(1)), 200); % time in seconds
Az_grid = linspace(min(Az_deg), max(Az_deg), 200);                  % azimuth degrees
[TIME_GRID, Az_GRID] = meshgrid(time_grid, Az_grid);

% Interpolate elevation data onto this grid
El_grid = griddata(seconds(time_vec - time_vec(1)), Az_deg, El_deg, TIME_GRID, Az_GRID, 'linear');

%% Plotting ISS Trajectory

% Create a figure
figure;

% Plot the ISS trajectory in 3D
plot3(r_ISS_ECI(:,1), r_ISS_ECI(:,2), r_ISS_ECI(:,3), 'LineWidth', 1.5);

% Set axis labels
xlabel('X [km]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
ylabel('Y [km]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
zlabel('Z [km]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');

% Add a title
title('Trajectory of the ISS in the ECI Frame', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');

% Add a grid
grid on;

% Make the plot square
axis equal;

% Set view for better visualization
view(3);

% Plot the Earth as a sphere at the origin (assuming Earth radius ~6371 km)
[XS, YS, ZS] = sphere(50);  % 50 is the resolution of the sphere

% Load the PNG image of Earth
earth_image = imread('C:\Users\scott\Documents\Folder\ME-594\Midterm\plots\earth_image.jpg');

% Plot Earth sphere with the texture
hold on;
surf(R_E*XS, R_E*YS, R_E*ZS, 'FaceColor', 'texturemap', 'CData', flipud(earth_image), 'EdgeColor', 'none');

% Plot the yellow star for r_0_s_ECI
h_star = plot3(x, y, z, 'y*', 'MarkerSize', 15, 'LineWidth', 3);

hold off;

% Adjust the axis limits to fit both the Earth and the ISS trajectory
axis([-1.1*R_E 1.1*R_E -1.1*R_E 1.1*R_E -1.1*R_E 1.1*R_E]);

% Ensure that the axis labels and the plot have consistent scaling
set(gca, 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');

% Add a legend for the sensor location (yellow star)
legend(h_star, 'Sensor Location: Albuquerque, New Mexico', 'FontSize', 14, 'FontWeight', 'bold');

%% Figure 2: Right Ascension vs Declination and Azimuth vs Elevation
figure;

% Subplot 1: Right Ascension vs Declination
subplot(1,2,1); % 1 row, 2 columns, plot 1
plot(RA_deg, Dec_deg, 'k.', 'LineWidth', 2)
xlabel('Right Ascension [deg]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold')
ylabel('Declination [deg]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold')
title('Right Ascension vs Declination', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold')
grid on;
axis square;

% Subplot 2: Azimuth vs Elevation
subplot(1,2,2); % 1 row, 2 columns, plot 2
plot(Az_deg, El_deg, 'k.', 'LineWidth', 2)
xlabel('Azimuth [deg]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold')
ylabel('Elevation [deg]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold')
title('Azimuth vs Elevation', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold')
grid on;
axis square;

%% Figure 3: Contour plot of Time vs Azimuth vs Elevation
figure;

contourf(TIME_GRID, Az_GRID, El_grid, 20, 'LineColor', 'none')
colorbar;
xlabel('Time', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold')
ylabel('Azimuth [deg]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold')
title('Elevation Contour Plot', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold')
grid on;
axis square;

% Adjust x-axis to display datetime labels
xticks = linspace(0, seconds(time_vec(end) - time_vec(1)), 8); % 8 ticks
xticklabels = datestr(date_0 + seconds(xticks), 'mmm dd HH:MM');
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold')

% Plot the point of maximum elevation as a red star
hold on;
az_max = rad2deg(AzEl_meas(max_el_idx,1));         % Azimuth in degrees
el_max_deg = rad2deg(max_el);                      % Max elevation in degrees
plot(t_max_el, az_max, 'r*', 'MarkerSize', 20, 'LineWidth', 2);  % Red star

% Display annotation text above the point with elevation value and date
date_at_max_el = datestr(date_0 + seconds(t_max_el), 'mmm dd, yyyy HH:MM:SS');  % Convert time to date string
text(t_max_el - 200000, az_max + 5, ...
    sprintf('Max Elevation\n%.2f°\nDate: %s', el_max_deg, date_at_max_el), ...
    'Color', 'black', 'FontSize', fs, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'FontName', 'Times New Roman');
hold off;