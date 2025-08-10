function az_el = compute_AzEl(X0, t0, tf, r_s_TH, L, lon, date)
% DESCRIPTION
%   Computes the azimuth and elevation of the International Space Station (ISS) 
%   relative to a specified ground station.
%
% INPUTS   Size     Type      Description                             Units
%   X0      (6,1)    (Double)   Initial state vector: [position; velocity]   [DU]
%   t0      (1,1)    (Double)   Initial time                               [TU]
%   tf      (1,1)    (Double)   Final time                                 [TU]
%   r_s_TH  (3,1)    (Double)   Position vector of the ground station       [DU]
%   L       (1,1)    (Double)   Latitude of the ground station              [deg]
%   lon     (1,1)    (Double)   Longitude of the ground station             [deg]
%   date    (1,1)    (Double)   Date for the calculation (Julian date)      [unitless]
%
% OUTPUTS  Size     Type      Description                             Units
%   az_el   (2,1)    (Double)   Azimuth and elevation of the ISS             [rad]
%         (1,1)     (Double)   Azimuth                                  [rad]
%         (2,1)     (Double)   Elevation                                [rad]
%
% FUNCTION
%   This function computes the azimuth and elevation of the International Space Station (ISS)
%   relative to a ground station. The calculation is done by propagating the ISS position using 
%   the two-body orbital model and then transforming the position into a topocentric coordinate 
%   system at the specified location of the ground station.
%
    [~, X_prop] = ode45(@two_body_ode, [t0 tf], X0);
    r_ISS = X_prop(end, 1:3)';

    [~, theta] = geocentric_to_ECI(r_s_TH, L, lon, date);
    Q_TH_ECI = rot_ijk_sez(L, theta);

    r_TH = Q_TH_ECI * r_ISS - r_s_TH;
    u_TH = r_TH / norm(r_TH);

    el = asin(-u_TH(3));
    az = mod(atan2(u_TH(1), u_TH(2)), 2*pi);

    if az < 0
        az = az + 2*pi;
    end

    az_el = [az; el];
end