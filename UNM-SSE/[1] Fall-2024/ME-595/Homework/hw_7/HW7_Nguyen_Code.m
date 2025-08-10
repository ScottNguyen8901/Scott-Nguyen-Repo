clear
clc
close all

% Define constants and parameters
mu = 1.3271244e11;  % Gravitational parameter of the Sun (km^3/s^2)
mu_p1 = 3.986e5;     % Gravitational parameter of Earth   (km^3/s^2)
mu_p2 = 4.2828e4;    % Gravitational parameter of Mars    (km^3/s^2)

R_p1 = 1.4959787e8;  % Distance from the Sun to Earth (km)
R_p2 = 2.279391e8;    % Distance from the Sun to Mars  (km)

rad_p1 = 6378;  % Radius of Earth (km)
rad_p2 = 3400;  % Radius of Mars  (km)

alt_p1 = 300;  % Altitude of parking orbit around Earth (km)
alt_p2 = 200;  % Altitude of parking orbit around Mars  (km)

% Calculate the radius of the parking orbits around Earth and Mars
r_p1 = rad_p1 + alt_p1;  % Radius of the parking orbit around Earth (km)
r_p2 = rad_p2 + alt_p2;    % Radius of the parking orbit around Mars  (km)

% Calculate the transfer orbit parameters
x = R_p2 / R_p1;  % Ratio of distances from the Sun to Mars and Earth
nu = linspace(pi/2, pi);  % Angle array from pi/2 to pi
% Calculate eccentricity for the transfer orbits
e = (x - 1) ./ (1 - x * cos(nu));
% Calculate semi-major axes for the transfer orbits
a = R_p1 * ((1 - x * cos(nu)) ./ (2 - x * (1 + cos(nu))));

% Calculate time of flight for each transfer orbit
tof = tof_ta(a, e, zeros(1, length(nu)), nu, mu);

% Perform patch conic transfer calculations
[dv_p1, dv_p2, dv_t, v_inf_p1, v_inf_p2, dv_t_hel] =...
    patch_conic_transfer(mu, mu_p1, mu_p2, R_p1, R_p2, r_p1, r_p2, a, e);

% Plot total Delta-V for the mission
figure;
plot(rad2deg(nu), dv_t, 'DisplayName', 'Total Delta-V');
hold on
plot(rad2deg(nu), dv_t_hel, 'DisplayName', 'Heliocentric Delta-V');
xlabel('Transfer Angle (degrees)');  % Label for x-axis
ylabel('Delta-V (km/s)');  % Label for y-axis
title('Delta-V vs. Transfer Angle');  % Title for the plot
legend show;  % Show legend for the plot

% Plot time of flight for the transfer orbits
figure;
plot(rad2deg(nu), tof / 86400);  % Convert time of flight from seconds to days
xlabel('Transfer Angle (degrees)');  % Label for x-axis
ylabel('Time of Flight (days)');  % Label for y-axis
title('Time of Flight vs. Transfer Angle');  % Title for the plot
