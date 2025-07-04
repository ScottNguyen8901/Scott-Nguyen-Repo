clear
close all
clc

% Populating constants
constants

% Define the orbital elements in km, radians, seconds
koe.a = 6785.6;              % Semi-major axis       [km]
koe.e = 0.01;                % Eccentricity          []
koe.i = deg2rad(35);         % Inclination           [rad]
koe.W = deg2rad(13.9687);    % RAAN                  [rad]
koe.w = deg2rad(82.7408);    % Argument of Periapsis [rad]
koe.f = deg2rad(281.9577);   % True anomaly          [rad]

% Convert orbital elements to state vector
x_0_ISS_ECI = koe_to_rv(koe, mu_E); % [km; km/s]

% Time span for integration
P = 2 * pi * sqrt(koe.a^3 / mu_E); % [s]
t_span = [0 20*P];                 % [s]

% Load atmosphere data
atmo_data = readmatrix('atmo_table.csv', 'NumHeaderLines', 1);

% Precompute atmospheric density values for a range of altitudes
altitudes = linspace(0, 1000, 1001);
rho_table = zeros(length(altitudes), 2);

for i = 1:length(altitudes)
    rho_table(i, 1) = altitudes(i);
    rho_table(i, 2) = comp_rho(altitudes(i), atmo_data);
end

%% Integrate the ODEs

options = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

[t, x_ISS_ECI]           =...
    ode45(@(t, y) two_body_ode(t, y), t_span, x_0_ISS_ECI, options);
[t_J2, x_ISS_ECI_J2]     =...
    ode45(@(t, y) two_body_ode_J2(t, y), t_span, x_0_ISS_ECI, options);
[t_drag, x_ISS_ECI_drag] =...
    ode45(@(t, y) two_body_ode_drag(t, y, rho_table), t_span, x_0_ISS_ECI, options);

%% Converting to Orbital Elements

oe      = rv_to_koe(x_ISS_ECI, mu_E);
oe_J2   = rv_to_koe(x_ISS_ECI_J2, mu_E);
oe_drag = rv_to_koe(x_ISS_ECI_drag, mu_E);

%% Plotting ISS Trajectories (With and Without J2 Perturbation, and With Drag)

fs = 24;

% ---- Plot Without J2 Perturbation ----
figure;
plot3(x_ISS_ECI(:,1), x_ISS_ECI(:,2), x_ISS_ECI(:,3), 'b', 'LineWidth', 2);  % Blue color with LineWidth 2
hold on;

% Plot Earth
[XS, YS, ZS] = sphere(50);
earth_image = imread('C:\Users\scott\Documents\Folder\ME-594\earth_image.jpg');
surf(R_E*XS, R_E*YS, R_E*ZS, 'FaceColor', 'texturemap', 'CData', flipud(earth_image), 'EdgeColor', 'none');

xlabel('X [km]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
ylabel('Y [km]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
zlabel('Z [km]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
title('Two-Body Trajectory', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');

axis equal;
axis([-1.1*R_E 1.1*R_E -1.1*R_E 1.1*R_E -1.1*R_E 1.1*R_E]);
grid on;
view(3);

% Create shared Earth surface mesh
[XS, YS, ZS] = sphere(50);

% ---- Combined Figure: J2 (left) and Drag (right) ----
figure('Name', 'Orbits with J2 and Drag Perturbations', 'Position', [100, 100, 1200, 600]);

% --- Left subplot: J2 ---
subplot(1,2,1);
plot3(x_ISS_ECI_J2(:,1), x_ISS_ECI_J2(:,2), x_ISS_ECI_J2(:,3), 'b', 'LineWidth', 2);
hold on;
surf(R_E*XS, R_E*YS, R_E*ZS, 'FaceColor', 'texturemap', 'CData', flipud(earth_image), 'EdgeColor', 'none');
xlabel('X [km]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
ylabel('Y [km]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
zlabel('Z [km]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
title('J2 Perturbation', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
axis equal;
axis([-1.1*R_E 1.1*R_E -1.1*R_E 1.1*R_E -1.1*R_E 1.1*R_E]);
grid on;
view(3);

% --- Right subplot: Drag ---
subplot(1,2,2);
plot3(x_ISS_ECI_drag(:,1), x_ISS_ECI_drag(:,2), x_ISS_ECI_drag(:,3), 'b', 'LineWidth', 2);
hold on;
surf(R_E*XS, R_E*YS, R_E*ZS, 'FaceColor', 'texturemap', 'CData', flipud(earth_image), 'EdgeColor', 'none');
xlabel('X [km]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
ylabel('Y [km]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
zlabel('Z [km]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
title('Drag Perturbation', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
axis equal;
axis([-1.1*R_E 1.1*R_E -1.1*R_E 1.1*R_E -1.1*R_E 1.1*R_E]);
grid on;
view(3);

%% Plotting Orbital Elements

% ---- Plot Orbital Elements Over Time Two-Body ----
figure('Name', 'Orbital Elements Over Time Two-Body', 'Position', [100, 100, 800, 1000]);
titles = {'Semi-major Axis [km]', 'Eccentricity', 'Inclination [deg]', ...
          'RAAN [deg]', 'Argument of Periapsis [deg]', 'True Anomaly [deg]'};

% Convert orbital elements from radians to degrees for relevant titles
oe(:,3:6) = rad2deg(oe(:,3:6));  % Inclination, RAAN, Argument of Periapsis, True Anomaly

% Divide time by orbital period P to express in periods
t_periods = t / P;

% Plot each orbital element
for k = 1:6
    subplot(6,1,k);
    
    % Plot the variable orbital element in black
    h1 = plot(t_periods, oe(:,k), 'k', 'LineWidth', 2);
    hold on;

    % Plot the constant orbital element in red (except for True Anomaly)
    if k ~= 6  % Exclude True Anomaly from the constant plot
        switch k
            case 1
                h2 = plot(t_periods, repmat(koe.a, size(t_periods)), 'r--', 'LineWidth', 2);  % Semi-major axis
            case 2
                h2 = plot(t_periods, repmat(koe.e, size(t_periods)), 'r--', 'LineWidth', 2);  % Eccentricity
            case 3
                h2 = plot(t_periods, repmat(rad2deg(koe.i), size(t_periods)), 'r--', 'LineWidth', 2);  % Inclination
            case 4
                h2 = plot(t_periods, repmat(rad2deg(koe.W), size(t_periods)), 'r--', 'LineWidth', 2);  % RAAN
            case 5
                h2 = plot(t_periods, repmat(rad2deg(koe.w), size(t_periods)), 'r--', 'LineWidth', 2);  % Argument of Periapsis
        end
    end

    % Use Greek letters for orbital elements in ylabel
    switch k
        case 1
            ylabel('$a$ [km]', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
        case 2
            ylabel('$e$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
        case 3
            ylabel('$i$ [deg]', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
        case 4
            ylabel('$\Omega$ [deg]', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
        case 5
            ylabel('$\omega$ [deg]', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
        case 6
            ylabel('$\nu$ [deg]', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
    end
    grid on;
    if k == 6
        xlabel('Time [Periods]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
    end
end

% Add a legend for the plot, specifying handles for the black and red lines
legend([h1, h2], {'Calculated', 'Constant'}, 'Location', 'Best', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

sgtitle('Orbital Elements: Two Body', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');  % Big title
% ---- Plot Comparison: J2 vs Drag Orbital Elements ----
figure('Name', 'J2 vs Drag Orbital Elements Comparison', 'Position', [100, 100, 1200, 1000]);

% Convert orbital elements from radians to degrees for relevant titles
oe_J2(:,3:6) = rad2deg(oe_J2(:,3:6));
oe_drag(:,3:6) = rad2deg(oe_drag(:,3:6));

% Divide time by orbital period P to express in periods
t_J2_periods = t_J2 / P;
t_drag_periods = t_drag / P;

for k = 1:6
    % Left column: J2 perturbation
    subplot(6,2,2*k-1);
    plot(t_J2_periods, oe_J2(:,k), 'k', 'LineWidth', 2);  % 'k' for black color
    % Use Greek letters for orbital elements in ylabel
    switch k
        case 1
            ylabel('$a$ [km]', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
        case 2
            ylabel('$e$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
        case 3
            ylabel('$i$ [deg]', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
        case 4
            ylabel('$\Omega$ [deg]', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
        case 5
            ylabel('$\omega$ [deg]', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
        case 6
            ylabel('$\nu$ [deg]', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
    end
    grid on;
    if k == 6
        xlabel('Time [Periods]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
    end

    % Right column: Drag perturbation
    subplot(6,2,2*k);
    plot(t_drag_periods, oe_drag(:,k), 'k', 'LineWidth', 2);  % 'k' for black color
    grid on;
    if k == 6
        xlabel('Time [Periods]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
    end
end
sgtitle('Orbital Elements: J2 vs Drag', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');  % Big title