clear
close all
clc

% Populating constants
constants

% ==========================
% Visibility criteria (NEW)
% ==========================
el_min      = deg2rad(10);     % minimum elevation for "viewable" [rad]
sunElNight  = deg2rad(-18);    % astronomical night threshold [rad] (-12 nautical, -6 civil)

R_Sun = 696340;               % Sun radius [km] (for shadow cone)
% ==========================

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
L   = deg2rad(35.0844);   % [rad]
lon = deg2rad(-106.6504); % [rad]

% Converting orbital elements to state vector
x_0_ISS_ECI = koe_to_rv(koe, mu_E); % [km; km/s]

% Specifying time since epoch to numerically integrate
tf = seconds(date_f - date_0);

% ==========================
% Increase fidelity (NEW)
% ==========================
% More samples => better chance of catching max elevation under constraints.
% Keep it reasonable so runtime doesn't explode.
numPts = max(5000, min(50000, round(tf/10)));  % ~1 sample per 10 seconds (clamped)
t_span = linspace(0, tf, numPts);
% ==========================

% Integrate the ODEs
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[t, x_ISS_ECI] = ode45(@(t, y) two_body_ode(t, y), t_span, x_0_ISS_ECI, options);

r_ISS_ECI = x_ISS_ECI(:,1:3); %Extracting ISS position

N = length(t);             % Number of measurements
LOS_meas   = zeros(N,2);   % Allocating matrix for RA/Dec measurements
AzEl_meas  = zeros(N,2);   % Allocating matrix for Az/El measurements

% ==========================
% Visibility bookkeeping (NEW)
% ==========================
sunEl_obs   = zeros(N,1);  % Sun elevation at observer [rad]
isEclipsed  = false(N,1);  % ISS in Earth's umbra? (conical shadow)
isNight     = false(N,1);  % observer nighttime?
isVisible   = false(N,1);  % final "viewable" mask
phaseAngle  = zeros(N,1);  % optional: Sun-ISS-Observer phase angle [rad]
% ==========================

for i = 1:N

    date = date_0 + seconds(t(i));                     % Current date
    theta_g = GMST(date);                              % GMST
    theta = lon + theta_g;

    Q_TH_ECI = rot_ijk_sez(L, theta);                  % ECI->TH (SEZ-ish)

    % --------------------------
    % ISS Az/El
    % --------------------------
    r_ISS_TH = Q_TH_ECI * r_ISS_ECI(i,:)';             % ISS position TH
    r_s_ISS_TH = r_ISS_TH - r_s_TH;                    % Sens to ISS in TH
    u_s_ISS_TH = r_s_ISS_TH / norm(r_s_ISS_TH);        % unit vec in TH

    el = asin(-u_s_ISS_TH(3));                         % Elevation: arcsin(-z)
    az = atan2(u_s_ISS_TH(1), u_s_ISS_TH(2));          % Azimuth: atan2(x, y)
    if az < 0
        az = az + 2*pi;
    end
    AzEl_meas(i,:) = [az, el];

    % --------------------------
    % RA/Dec
    % --------------------------
    r_s_ECI = Q_TH_ECI' * r_s_TH;                      % sensor position in ECI
    r_s_ISS_ECI = r_ISS_ECI(i,:)' - r_s_ECI;           % sensor -> ISS
    u_s_ISS_ECI = r_s_ISS_ECI / norm(r_s_ISS_ECI);

    RA  = atan2(u_s_ISS_ECI(2), u_s_ISS_ECI(1));
    if RA < 0
        RA = RA + 2*pi;
    end
    Dec = asin(u_s_ISS_ECI(3));
    LOS_meas(i,:) = [RA, Dec];

    % ==========================================================
    % Sun geometry + visibility conditions (NEW, higher fidelity)
    % ==========================================================

    % Sun position in ECI (Earth-centered) [km]
    r_Sun_ECI = sun_eci_approx(date);
    r_ES = r_Sun_ECI;                         % Earth -> Sun
    d_ES = norm(r_ES);
    u_ES = r_ES / d_ES;                       % unit Earth->Sun

    % 1) Observer nighttime: Sun elevation at observer
    r_Sun_TH = Q_TH_ECI * r_Sun_ECI;          % Sun position in TH
    r_s_Sun_TH = r_Sun_TH - r_s_TH;           % sensor -> Sun (TH)
    u_s_Sun_TH = r_s_Sun_TH / norm(r_s_Sun_TH);

    sunEl_obs(i) = asin(-u_s_Sun_TH(3));      % Sun elevation in same convention
    isNight(i)   = (sunEl_obs(i) < sunElNight);

    % 2) ISS sunlit? Conical umbra test (Earth shadow cone)
    r_IS = r_ISS_ECI(i,:)';                   % Earth -> ISS

    % Components relative to Sun-line
    proj = dot(r_IS, u_ES);                   % along +Sun direction
    x = -proj;                                % distance "behind Earth" along anti-sun axis

    if x <= 0
        % Not behind Earth relative to Sun => cannot be in umbra
        isEclipsed(i) = false;
    else
        % Perpendicular distance from Sun-line axis
        d_perp = norm(r_IS - proj*u_ES);

        % Umbra cone radius at distance x behind Earth:
        % R_umbra(x) = R_E - x*(R_Sun - R_E)/d_ES
        R_umbra = R_E - x*(R_Sun - R_E)/d_ES;

        % If cone radius is still positive and sat lies within it => eclipsed
        isEclipsed(i) = (R_umbra > 0) && (d_perp < R_umbra);
    end

    % 3) Optional: phase angle (Sun-ISS-Observer angle at ISS)
    u_ISS_to_Sun = (r_Sun_ECI - r_IS); u_ISS_to_Sun = u_ISS_to_Sun / norm(u_ISS_to_Sun);
    u_ISS_to_Obs = (r_s_ECI   - r_IS); u_ISS_to_Obs = u_ISS_to_Obs / norm(u_ISS_to_Obs);
    phaseAngle(i) = acos( max(-1, min(1, dot(u_ISS_to_Sun, u_ISS_to_Obs))) );

    % Final visibility mask (AND of constraints)
    isVisible(i) = (el > el_min) && isNight(i) && (~isEclipsed(i));
    % ==========================================================
end

% ==========================================================
% ORIGINAL: Find the maximum elevation and its corresponding time
% (KEEPING YOUR ORIGINAL PRINTING)
% ==========================================================
[max_el, max_el_idx] = max(AzEl_meas(:,2));  % max elevation and its index
t_max_el = t(max_el_idx);                    % time at max elevation

% Display result
fprintf('Maximum Elevation: %.6f radians (%.3f deg)\n', max_el, rad2deg(max_el));
fprintf('Time since Epoch at Max Elevation: %.3f seconds\n', t_max_el);

% Extract the corresponding state from x_ISS_ECI
state_at_dt = x_ISS_ECI(max_el_idx, :);

% Extract position and velocity components
position = state_at_dt(1:3);
velocity = state_at_dt(4:6);

% Extract azimuth and elevation at dt
RA_Dec_at_dt = LOS_meas(max_el_idx, :);  % [Azimuth, Elevation] at dt (your original comment)
RA = RA_Dec_at_dt(1);                    % Right Ascension in radians
Dec = RA_Dec_at_dt(2);                   % Declination in radians

% Calculate the corresponding date for dt
date_at_dt = date_0 + seconds(t_max_el);

% Display the state vector in a neat format
disp('State at dt:');
disp(' ');  % Adding a space for better readability

% Display position
fprintf('Position (x, y, z) [km]:\n');
fprintf('x = %.4f km\n', position(1));
fprintf('y = %.4f km\n', position(2));
fprintf('z = %.4f km\n', position(3));

disp(' ');

% Display velocity
fprintf('Velocity (vx, vy, vz) [km/s]:\n');
fprintf('vx = %.4f km/s\n', velocity(1));
fprintf('vy = %.4f km/s\n', velocity(2));
fprintf('vz = %.4f km/s\n', velocity(3));

disp(' ');

% Display azimuth and elevation
fprintf('Right Ascension and Declination [deg]:\n');
fprintf('Right Ascension = %.4f rad\n', RA);
fprintf('Declination = %.4f rad\n', Dec);

disp(' ');

% Display the corresponding date
fprintf('Corresponding date: %s\n', datestr(date_at_dt));

% ==========================================================
% NEW: Best viewing time = max elevation subject to visibility constraints
% ==========================================================
if any(isVisible)
    vis_idx = find(isVisible);
    [best_el, kbest] = max(AzEl_meas(vis_idx,2));
    best_idx = vis_idx(kbest);

    t_best = t(best_idx);
    date_best = date_0 + seconds(t_best);

    fprintf('\n=== BEST VIEWING TIME (max elevation with all criteria true) ===\n');
    fprintf('Elevation: %.6f radians (%.3f deg)\n', best_el, rad2deg(best_el));
    fprintf('Azimuth:   %.6f radians (%.3f deg)\n', AzEl_meas(best_idx,1), rad2deg(AzEl_meas(best_idx,1)));
    fprintf('Time since Epoch: %.3f seconds\n', t_best);
    fprintf('Corresponding date: %s\n', datestr(date_best));

    fprintf('Sun elevation at observer: %.3f deg\n', rad2deg(sunEl_obs(best_idx)));
    fprintf('ISS eclipsed? %d\n', isEclipsed(best_idx));
else
    fprintf('\nNo times satisfy visibility criteria in this window.\n');
end

% ==========================
% (Optional) List visible time windows (kept)
% ==========================
time_vec = date_0 + seconds(t);
vis = isVisible(:);
edges = diff([false; vis; false]);
startIdx = find(edges == 1);
endIdx   = find(edges == -1) - 1;

if ~isempty(startIdx)
    fprintf('\nVisible windows (el>%.1f deg, SunEl<%.1f deg, not eclipsed):\n', ...
        rad2deg(el_min), rad2deg(sunElNight));
    for k = 1:numel(startIdx)
        fprintf('  %s  ->  %s\n', datestr(time_vec(startIdx(k))), datestr(time_vec(endIdx(k))));
    end
end

%% Plotting

% Set the desired font size
fs = 24;

% Convert Right Ascension and Declination to degrees
RA_deg = rad2deg(LOS_meas(:,1));
Dec_deg = rad2deg(LOS_meas(:,2));

% Convert Azimuth and Elevation to degrees
Az_deg = rad2deg(AzEl_meas(:,1));
El_deg = rad2deg(AzEl_meas(:,2));

% Create grids for interpolation
time_grid = linspace(0, seconds(time_vec(end) - time_vec(1)), 200);
Az_grid = linspace(min(Az_deg), max(Az_deg), 200);
[TIME_GRID, Az_GRID] = meshgrid(time_grid, Az_grid);

% Interpolate elevation data onto this grid
El_grid = griddata(seconds(time_vec - time_vec(1)), Az_deg, El_deg, TIME_GRID, Az_GRID, 'linear');

%% Plotting ISS Trajectory
x = R_E * cos(L) * cos(lon);
y = R_E * cos(L) * sin(lon);
z = R_E * sin(L);

figure;
plot3(r_ISS_ECI(:,1), r_ISS_ECI(:,2), r_ISS_ECI(:,3), 'LineWidth', 2);
xlabel('X [km]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
ylabel('Y [km]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
zlabel('Z [km]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
title('Trajectory of the ISS in the ECI Frame', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
grid on; axis equal; view(3);

[XS, YS, ZS] = sphere(50);
earth_image = imread("C:\Users\scott\Documents\Folder\UNM-SSE\2-Spring-2025\ME-594\earth_image.jpg");

hold on;
surf(R_E*XS, R_E*YS, R_E*ZS, 'FaceColor', 'texturemap', 'CData', flipud(earth_image), 'EdgeColor', 'none');
h_star = plot3(x, y, z, 'y*', 'MarkerSize', 15, 'LineWidth', 3);
hold off;

axis([-1.1*R_E 1.1*R_E -1.1*R_E 1.1*R_E -1.1*R_E 1.1*R_E]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold');
legend(h_star, 'Sensor Location: Albuquerque, New Mexico', 'FontSize', 14, 'FontWeight', 'bold');

%% Figure 2: Right Ascension vs Declination and Azimuth vs Elevation (NO LEGENDS)
figure;

subplot(1,2,1);
plot(RA_deg, Dec_deg, 'k.', 'LineWidth', 2); hold on;
plot(RA_deg(isVisible), Dec_deg(isVisible), 'r.', 'LineWidth', 2); hold off;
xlabel('Right Ascension [deg]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold')
ylabel('Declination [deg]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold')
title('Right Ascension vs Declination', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold')
grid on; axis square;
text(0.02,0.95,'Red = Visible','Units','normalized','FontName','Times New Roman', ...
     'FontSize',14,'FontWeight','bold');

subplot(1,2,2);
plot(Az_deg, El_deg, 'k.', 'LineWidth', 2); hold on;
plot(Az_deg(isVisible), El_deg(isVisible), 'r.', 'LineWidth', 2); hold off;
xlabel('Azimuth [deg]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold')
ylabel('Elevation [deg]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold')
title('Azimuth vs Elevation', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold')
grid on; axis square;
text(0.02,0.95,'Red = Visible','Units','normalized','FontName','Times New Roman', ...
     'FontSize',14,'FontWeight','bold');

%% Figure 3: Contour plot + mark BEST VIEWABLE time only
figure;
contourf(TIME_GRID, Az_GRID, El_grid, 20, 'LineColor', 'none')
colorbar;
xlabel('Time', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold')
ylabel('Azimuth [deg]', 'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold')
title('Elevation Contour Plot (Best Viewable Time)', ...
      'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold')
grid on; axis square;

% Time axis formatting
xticks = linspace(0, seconds(time_vec(end) - time_vec(1)), 8);
xticklabels = datestr(date_0 + seconds(xticks), 'mmm dd HH:MM');
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, ...
         'FontName', 'Times New Roman', 'FontSize', fs, 'FontWeight', 'bold')

hold on;

% ==========================================================
% BEST VIEWABLE TIME (passes ALL criteria)
% ==========================================================
if any(isVisible)
    vis_idx = find(isVisible);
    [~, kbest] = max(AzEl_meas(vis_idx,2));   % max elevation among visible
    best_idx = vis_idx(kbest);

    % Azimuth at best time
    az_best = rad2deg(AzEl_meas(best_idx,1));

    % Plot red star
    plot(t(best_idx), az_best, 'r*', 'MarkerSize', 22, 'LineWidth', 2);

    % Local time label
    best_time_local = date_0 + seconds(t(best_idx));
    label_str = datestr(best_time_local, 'mmm dd, HH:MM');

    text(t(best_idx) - 5e4, az_best - 20, ...
        sprintf('Best Viewable\n%s', label_str), ...
        'FontName', 'Times New Roman', ...
        'FontSize', 16, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center');
end

hold off;
