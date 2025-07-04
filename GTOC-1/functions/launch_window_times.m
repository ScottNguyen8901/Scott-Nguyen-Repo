function [LW_AN, LW_DN] = launch_window_times(tle_epoch, i, W, w, M, e, n, lat, lon)
    %
    % DESCRIPTION
    %   Calculate launch window times at the ascending and descending nodes based on the
    %   orbital elements and launch site location. This function computes the times when 
    %   the satellite passes over the launch site, providing two possible launch windows.
    %
    % INPUTS         size         Type    Description                        Units
    %   tle_epoch     (1,1)       Double  TLE epoch time (decimal days)      [days]
    %   i             (1,1)       Double  Orbital inclination                [deg]
    %   W             (1,1)       Double  Right Ascension of Ascending Node  [deg]
    %   w             (1,1)       Double  Argument of perigee                [deg]
    %   M             (1,1)       Double  Mean anomaly                       [deg]
    %   e             (1,1)       Double  Orbital eccentricity               []
    %   n             (1,1)       Double  Mean motion (revolutions/day)      [rev/day]
    %   lat           (1,1)       Double  Latitude of the launch site        [deg]
    %   lon           (1,1)       Double  Longitude of the launch site       [deg]
    %
    % OUTPUTS        size         Type    Description                        Units
    %   LW_AN         (1,1)       Datetime Launch window at ascending node   [datetime]
    %   LW_DN         (1,1)       Datetime Launch window at descending node  [datetime]
    %
    % NOTES
    %   This function assumes a circular Earth model and uses basic Keplerian orbital mechanics
    %   to compute launch windows. The function does not account for perturbations or 
    %   other dynamic effects that could affect the launch window timings.
    %
    % FUNCTION

    % Constants
    mu = 398600;            % Earth's gravitational parameter in km^3/s^2
    w_earth = 7.2921159e-5; % Earth angular velocity in rad/s

    % Convert input angles from degrees to radians
    L_s = deg2rad(lat);      % Launch site latitude
    lambda = deg2rad(lon);   % Launch site longitude (east +, west -)

    % Convert TLE epoch to Julian date
    tle_epoch_date_vec = tle_epoch_2_datetime(tle_epoch);
    tle_epoch_JD = juliandate(tle_epoch_date_vec);

    % Convert orbital elements from degrees to radians
    i = deg2rad(i);  
    W = deg2rad(W);          
    w = deg2rad(w);   
    M = deg2rad(M);  

    % Orbital parameters
    e = e;
    n = n * (2 * pi / 86400);  % Mean motion in rad/s
    a = (mu / n^2)^(1 / 3);  % Semi-major axis

    % Launch direction angles
    alpha = i;
    gamma = asin(cos(alpha) / cos(L_s));
    delta = acos(cos(gamma) / sin(alpha));

    % Ascending and descending node opportunities
    LWST_AN = W + delta;
    beta_AN = gamma;

    LWST_DN = W + pi - delta;
    beta_DN = pi - gamma;

    % Calculate sidereal time at epoch
    [theta_g, ~] = siderealTime(tle_epoch_JD);

    % Local sidereal time
    theta = deg2rad(theta_g) + lambda;

    % Wrap angles between 0 and 2*pi
    LWST_AN = mod(LWST_AN, 2*pi);
    LWST_DN = mod(LWST_DN, 2*pi);
    theta = mod(theta, 2*pi);

    % Calculate wait times at nodes
    delta_t_AN = (LWST_AN - theta) / w_earth;
    delta_t_DN = (LWST_DN - theta) / w_earth;

    % Calculate launch windows at nodes
    LW_AN = datetime(tle_epoch_JD + delta_t_AN / 86400, 'ConvertFrom', 'juliandate');
    LW_DN = datetime(tle_epoch_JD + delta_t_DN / 86400, 'ConvertFrom', 'juliandate');
end
