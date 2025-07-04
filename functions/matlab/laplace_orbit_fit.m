function x = laplace_orbit_fit(r_s_TH, lat, lon, T, RA_DEC)
%
% DESCRIPTION
%   Estimates the ECI (Earth-Centered Inertial) state vector using 
%   Laplace's method based on Right Ascension (RA) and Declination (DEC) 
%   observations from a ground station.
%
% INPUTS    size     Type     Description           Units
%   lat     (1,1)    Double   Latitude              [rad]
%   lon     (1,1)    Double   Longitude             [rad]
%   T       (1,1)    Double   Measurement Times     [TU]
%   RA_DEC  (3,2)    Double   RA and DEC            [rad]
%
% OUTPUTS   size     Type     Description                    Units
%   x       (6,1)    Double   Estimated ECI state at T_2     [DU, DU/TU]
%
% NOTES
%   The input RA/DEC are assumed to be in radians and the function 
%   calculates the ECI state vector at the middle of the observation times.
%
% FUNCTION

    constants; % Loading in constants

    % Reformat RA_DEC to match old format
    RA  = RA_DEC(:,1); % RA1, RA2, RA3
    DEC = RA_DEC(:,2); % DEC1, DEC2, DEC3

    % Build line-of-sight unit vectors
    L = [cos(RA) .* cos(DEC), ...
         sin(RA) .* cos(DEC), ...
         sin(DEC)]';
    
    % Converting times to Julian Date
    t_1 = juliandate(T(1));
    t_2 = juliandate(T(2));
    t_3 = juliandate(T(3));
    t   = juliandate(T(2));

    % Calculating first and second derivatives
    L_d = (2*t - t_2 - t_3)/((t_1 - t_2)*(t_1 - t_3)) * L(:,1) +...
          (2*t - t_1 - t_3)/((t_2 - t_1)*(t_2 - t_3)) * L(:,2) +...
          (2*t - t_1 - t_2)/((t_3 - t_1)*(t_3 - t_2)) * L(:,3);

    L_dd = 2/((t_1 - t_2)*(t_1 - t_3)) * L(:,1) +...
          2/((t_2 - t_1)*(t_2 - t_3)) * L(:,2) +...
          2/((t_3 - t_1)*(t_3 - t_2)) * L(:,3);
    
    % Scaling by time
    L_d = L_d / sec_day;
    L_dd = L_dd / sec_day^2;

    % Get ground station ECI position and velocity at T(2)
    [x_s_2_ECI, ~] = geocentric_to_ECI(r_s_TH, lat, lon, T(2));
    
    % Calculate terms needed for Laplace method
    R = norm(x_s_2_ECI(1:3));
    N = dot(L(:,2), x_s_2_ECI(1:3));
    D = 2 * det([L(:,2), L_d, L_dd]);
    D1 = det([L(:,2), L_d, x_s_2_ECI(7:9)]);
    D2 = det([L(:,2), L_d, x_s_2_ECI(1:3)]);
    D3 = det([L(:,2), x_s_2_ECI(7:9), L_dd]);
    D4 = det([L(:,2), x_s_2_ECI(1:3), L_dd]);

    % Coefficients for 8th degree polynomial
    c8 = D^2;
    c6 = (4*N*D - 4*D1)*D1 - D^2 * R^2;
    c3 = 4*mu_E*D2*(N*D - 2*D1);
    c0 = -4*mu_E^2*D2^2;

    % Solve for roots
    Z = roots([c8, 0, c6, 0, 0, c3, 0, 0, c0]);

    % Pick best positive real root
    r = Z(4);
    rrr = r^3;

    % Compute rho and rhodot
    rho = -2*(D1/D) - 2*(mu_E/rrr)*(D2/D);
    rhodot = -(D3/D) - (mu_E/rrr)*(D4/D);

    % Compute ECI position and velocity
    x = [rho * L(:,2) + x_s_2_ECI(1:3);
          rhodot * L(:,2) + rho * L_d + x_s_2_ECI(4:6)];
end