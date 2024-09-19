function [f_p, df_dp] = tof_gauss_lambert(p, r_1_vec, r_2_vec, dt, mu, orbit)
    %
    % DESCRIPTION
    %   Compute the time of flight between two position vectors using the Gauss-Lambert method.
    %   The function takes the variable 'p', position vectors at the initial and final points,
    %   time of flight, gravitational parameter, and the orbit type (short or long) as inputs to
    %   compute the residual for the time of flight.
    %
    % INPUTS         size         Type    Description                   Units
    %   p            (1,1)        Double  Iterative variable            []
    %   r_1_vec      (3,1)        Double  Initial position vector       [DU]
    %   r_2_vec      (3,1)        Double  Final position vector         [DU]
    %   dt           (1,1)        Double  Desired time of flight        [TU]
    %   mu           (1,1)        Double  Gravitational parameter       [DU^3/TU^2]
    %   orbit        (1,1)        Char    Orbit type ('short'/'long')   []
    %
    % OUTPUTS        size         Type    Description                   Units
    %   f_p          (1,1)        Double  Residual between desired TOF  [TU]
    %   df_dp        (1,1)        Double  Derivative of f_p w.r.t p     [DU/TU]
    %
    % NOTES
    %
    % FUNCTION

    r_1 = norm(r_1_vec);
    r_2 = norm(r_2_vec);
    
    switch orbit
        case 'short'
            delta_nu = acos(dot(r_1_vec, r_2_vec) / (r_1 * r_2));
        case 'long'
            delta_nu = 2*pi - acos(dot(r_1_vec, r_2_vec) / (r_1 * r_2));
    end
    
    k = r_1 * r_2 * (1 - cos(delta_nu));
    l = r_1 + r_2;
    m = r_1 * r_2 * (1 + cos(delta_nu));
    
    a = (m * k *p) / ((2 * m - l^2)*p^2 + 2 * k * l * p - k^2);
    
    f = 1 - (r_2 / p) * (1 - cos(delta_nu));
    delta_E = acos(1 - (r_1 / a) * (1 - f));
    g = r_1 * r_2 * sin(delta_nu) / sqrt(mu * p);
    
    dt_p = g + sqrt(a^3 / mu)*(delta_E - sin(delta_E));
    f_p = dt - dt_p;

    dt_dp = -(g/(2 * p)) - ...
            (3 / 2) * a * (dt_p - g) * ((k^2 + (2 * m - l^2) * p^2) / (m * k * p^2)) + ...
            sqrt(a^3 / mu) * ((2 * k * sin(delta_E)) / (mu * p * (k - l * p)));

    df_dp = -dt_dp;

end