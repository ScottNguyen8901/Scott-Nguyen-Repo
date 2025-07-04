clear
close all
clc

% Populating constants
constants

%% Problem 1

% Position location
lat     = deg2rad(35);   % [rad] 
lon     = deg2rad(-108); % [rad]
theta_g = deg2rad(5);    % [rad]
theta   = lon + theta_g; % [rad]

% Position of object in ECI
r_o_ECI = [5294.35;
           3707.14;
           2352.42]; % [km]

% Rotation matrix from ECI to SEZ
Q_ECI_TH = rot_eci_sez(lat, theta); % []

% Position of object in sensor TH frame
r_o_TH = Q_ECI_TH * r_o_ECI; % [km]

% Vector going from sensor to object in TH frame
r_so_TH = r_o_TH - r_s_TH; % [km]

% Printing
fprintf('Problem 1:\n\n');
disp('r_so_TH [km]:');
disp(r_so_TH);
fprintf('__________________________________________________________________\n\n');

%% Problem 2

% Right ascension and declination angles
RA  = deg2rad(17); % [rad]
Dec = deg2rad(62); % [rad]

% Unit vector to star in ECI
u_star_ECI = [cos(Dec)*cos(RA);
              cos(Dec)*sin(RA);
              sin(Dec)]; % [km]

% Unit vector to star in TH frame
u_star_TH = Q_ECI_TH * u_star_ECI; % [km]

% Printing
fprintf('Problem 2:\n\n');
disp('u_star_TH [km]:');
disp(u_star_TH);
