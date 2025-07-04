clear
close all
clc

% Populating constants
constants

%% Problem 1

% Camera specs
f     = 72;    % Focal length [mm]
cam_W = 1500;  % Camer width  [pixels]
cam_H = 2000;  % Camera heigt [pixels]
pix_W = 0.005; % Pixel width  [mm]
pix_H = pix_W; % Pixel height [mm]

% Position location
lat     = deg2rad(35);   % [rad] 
lon     = deg2rad(35);   % [rad]
theta_g = deg2rad(5);    % [rad]
theta   = lon + theta_g; % [rad]

% Position of objects in ECI
r_o_1_ECI = [5294.35;
             3707.14;
             2352.42]; % [km]
r_o_2_ECI = [5467.25;
             3489.87;
             2117.82]; % [km]

% Rotation matrix from ECI to SEZ
Q_ECI_TH = rot_eci_sez(lat, theta); % []

% Position of objects in sensor TH frame
r_o_1_TH = Q_ECI_TH * r_o_1_ECI; % [km]
r_o_2_TH = Q_ECI_TH * r_o_2_ECI; % [km]

% Vector going from sensor to objects in TH frame
r_so_1_TH = r_o_1_TH - r_s_TH; % [km]
r_so_2_TH = r_o_2_TH - r_s_TH; % [km]

% Calculate camera width and height FOV
[FOV_W, FOV_H] = calc_FOV(f, cam_W, cam_H, pix_W, pix_H); % [rad]

% Converting to degrees and finding if both objects can be in the FOV
FOV_W_deg = rad2deg(FOV_W); % [deg]
FOV_H_deg = rad2deg(FOV_H); % [deg]

alpha_deg = rad2deg(acos(dot(r_so_1_TH, r_so_2_TH) /...
    (norm(r_so_1_TH) * norm(r_so_2_TH)))); % [deg]

% Create and display the results as a table
T = table(FOV_W_deg, FOV_H_deg, alpha_deg, ...
    'VariableNames', {'FOV_Width_deg', 'FOV_Height_deg', 'Angle_Between_Objects_deg'});

% Printing
fprintf('Problem 1:\n\n');
disp(T);
fprintf('_____________________________________________________________________\n\n');

%% Problem 2

% Sensor position
lat = deg2rad(39);     % [rad]
lon = deg2rad(-104);   % [rad]
theta_g = 0;           % [rad]
theta = lon + theta_g; % [rad]

% Object range and azimuth and elevation
rho = 650.75;        % [km]
Az = deg2rad(25.12); % [rad]
El = deg2rad(17.29); % [rad]

% Vector from sensor to object in TH frame
r_so_TH = rho * [cos(El) * cos(Az);
                 cos(El) * sin(Az);
                 sin(El)];          % [km]

% Position of object in TH frame
r_o_TH = r_s_TH + r_so_TH; % [km]

% Rotation matrix from ECI to TH frame
Q_ECI_TH = rot_eci_sez(lat, theta); % []

% Position of object in ECI frame
r_o_ECI = Q_ECI_TH' * r_o_TH; % [km]

% Printing
fprintf('Problem 2:\n\n');
disp('r_o_ECI [km]:');
disp(r_o_ECI);