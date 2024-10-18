clear
clc
close all

% Define constants and parameters
mu_sun = 1.3271244e11;  % Gravitational parameter of the Sun (km^3/s^2)
mu_earth = 3.986e5;     % Gravitational parameter of Earth   (km^3/s^2)
mu_mars = 4.2828e4;     % Gravitational parameter of Mars    (km^3/s^2)

R_sun_earth = 1.4959787e8;  % Distance from the Sun to Earth (km)
R_sun_mars = 2.279391e8;    % Distance from the Sun to Mars  (km)

R_earth = 6378;  % Radius of Earth (km)
R_mars  = 3400;  % Radius of Mars  (km)

alt_earth = 300;  % Altitude of parking orbit around Earth (km)
alt_mars  = 200;  % Altitude of parking orbit around Mars  (km)

r_p_earth = R_earth + alt_earth;  % Radius of the parking orbit around Earth (km)
r_p_mars  = R_mars + alt_mars;    % Radius of the parking orbit around Mars  (km)

% Define eccentricity range and semi-latus rectum range (p = r(1 + e))
e_range = linspace(0.207, 0.523, 100);
p_range = R_sun_earth * (1 + e_range);  % Semi-latus rectum for each e
n = length(e_range);

% Matrices to store delta-v and time of flight
dv_total_pc = zeros(1, n);
dv_total_gt = zeros(1, n);
tof = zeros(1, n);

% True anomaly calculations
nu_1_range = real(acosd((p_range - R_sun_earth) ./ (e_range * R_sun_earth)));
nu_2_range = real(acosd((p_range - R_sun_mars) ./ (e_range * R_sun_mars)));

for i = 1:n
    % Extract eccentricity, semi-latus rectum and true anomaly
    e = e_range(i);
    p = p_range(i);
    nu_2 = nu_2_range(i);   
    
    % Call to functions to calculate delta-v
    [~, ~, dv_t_pc] = patch_conic_transfer(mu_sun, mu_earth, mu_mars, ...
        R_sun_earth, R_sun_mars, r_p_earth, r_p_mars, p, e);
    [~, ~, dv_t_gt] = general_transfer(R_sun_earth, R_sun_mars, p, e, mu_sun);
    
    % Storing delta-v and time of flight
    dv_total_pc(i) = dv_t_pc;
    dv_total_gt(i) = dv_t_gt;

    a_t = p / (1 - e^2);        % Semi-major axis of transfer orbit
    tof(i) = tof_ta(a_t, e, 0, deg2rad(nu_2), mu_sun);
end

% Converting time of flight to days
tof_days = tof / 86400;

% Plotting the results
figure;
plot(nu_2_range, dv_total_pc, 'b-', 'LineWidth', 2);
grid on;
axis square;
xlim([90 180]);
xlabel('\Delta\nu [degrees]', 'FontSize', 12);
ylabel('Total \Delta V [km/s]', 'FontSize', 12);
title('Patched Conics Transfer', 'FontSize', 14);
set(gca, 'FontSize', 12);

figure;
plot(nu_2_range, dv_total_gt, 'b-', 'LineWidth', 2);
grid on;
axis square;
xlim([90 180]);
xlabel('\Delta\nu [degrees]', 'FontSize', 12);
ylabel('Total \Delta V [km/s]', 'FontSize', 12);
title('Heliocentric Transfer', 'FontSize', 14);
set(gca, 'FontSize', 12);

figure;
plot(nu_2_range, tof_days, 'b-', 'LineWidth', 2);
grid on;
axis square;
xlim([90 180]);
xlabel('\Delta\nu [degrees]', 'FontSize', 12);
ylabel('Time of Flight [days]', 'FontSize', 12);
title('Time of Flight vs. True Anomaly', 'FontSize', 14);
set(gca, 'FontSize', 12);