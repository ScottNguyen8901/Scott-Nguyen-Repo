clear
close 
clc

% Given parameters
T_daylight = 68;  % minutes
T_eclipse = 22;  % minutes
P_daylight = 400;  % watts
P_eclipse = 100;  % watts
eta_daylight = 0.85;
eta_eclipse = 0.65;

% Power required from the solar array during daylight periods
P_solar_array = ((P_daylight * T_daylight / eta_daylight) + (P_eclipse * T_eclipse / eta_eclipse)) / T_daylight;

% Output the result with two decimal places
fprintf('Problem 1: Power required from the solar array during daylight periods: %.2f W\n', P_solar_array);

S = 1367;
eta = 1;
theta = 23.5;
I_d = 0.72;

% Beginning-of-life (BOL) power production capability required
P_BOL = P_BOL / (I_d * cosd(theta));

% Output the BOL power production capability
fprintf('Problem 2: Beginning-of-life power production capability required: %.2f W\n', P_BOL);

% BOL power per unit area
n = 0.24;  % Solar cell efficiency
Ip = 0.95;  % Packing factor
Itemp = 0.9;  % Temperature loss
Is = 1;  % No shadowing losses

P_BOL_per_m2 = solar_constant * n * Ip * Itemp * Is * cosd(incidence_angle);

% Output the BOL power per unit area
fprintf('Problem 3: Beginning-of-life power per unit area: %.2f W/m2\n', P_BOL_per_m2);

% Required solar array area
A_solar_array = P_BOL / P_BOL_per_m2;

% Output the required solar array area
fprintf('Problem 4: Required area of the solar array: %.2f m2\n', A_solar_array);

% Battery parameters
DOD = 0.20;  % Depth of discharge
eta_battery = 0.90;  % Transmission efficiency

% Energy required during eclipse
E_eclipse = P_eclipse * T_eclipse / eta_battery;

% Battery capacity required
battery_capacity = E_eclipse / DOD;

% Output the required battery capacity
fprintf('Problem 5: Required battery capacity: %.2f Wh\n', battery_capacity);

%%
clear
close
clc

% Given values
k = 0.005; % thermal conductivity in W/(m·K)
T_in = 20; % inner temperature in K
T_out = 300; % outer temperature in K
r_in = 0.125; % inner radius in meters
r_out = 0.145; % outer radius in meters

% Calculating the heat transfer rate (Q)
Q = (4 * pi * k * (T_out - T_in) * r_in * r_out) / log(r_out / r_in);

% Display the result
fprintf('Heat transfer rate: %.2f W\n', -Q);
