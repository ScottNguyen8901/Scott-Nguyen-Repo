function [v_1_vec, v_2_vec] = lambert_solver(root_solver, r_1_vec, r_2_vec, dt, mu, orbit)
%
% DESCRIPTION
%   Solve Lambert's problem for a two-body orbital transfer between two position vectors
%   and a given time of flight. The function uses the specified root-solving method 
%   (bisection or linear) to find the semilatus rectum (p), semi-major axis (a), 
%   Lagrange coefficients (f, g, f_dot, g_dot), and computes the initial and final
%   velocity vectors (v_1_vec and v_2_vec).
%
% INPUTS         size         Type    Description                               Units
%   root_solver  (1,1)        Char    Root solving method: 'bisection' or       [n/a]
%                                     'linear'
%   r_1_vec      (3,1)        Double  Initial position vector                   [km]
%   r_2_vec      (3,1)        Double  Final position vector                     [km]
%   dt           (1,1)        Double  Time of flight                            [s]
%   mu           (1,1)        Double  Gravitational parameter                   [km^3/s^2]
%   orbit        (1,1)        Char    Orbit type: 'short' or 'long' transfer    [n/a]
%
% OUTPUTS        size         Type    Description                               Units
%   v_1_vec      (3,1)        Double  Initial velocity vector                   [km/s]
%   v_2_vec      (3,1)        Double  Final velocity vector                     [km/s]
%
% NOTES
%   This function solves Lambert's problem using the f and g method. It computes 
%   the initial and final velocity vectors given the position vectors and time of flight.
%   The function supports two methods for root solving: linear and bisection. 
%

    % Normalize position vectors
    r_1 = norm(r_1_vec);
    r_2 = norm(r_2_vec);
    
    % Determine the change in true anomaly based on orbit type
    switch orbit
        case 'short'
            delta_nu = acos(dot(r_1_vec, r_2_vec) / (r_1 * r_2));
        case 'long'
            delta_nu = 2*pi - acos(dot(r_1_vec, r_2_vec) / (r_1 * r_2));
    end

    % Compute constants for solving
    k = r_1 * r_2 * (1 - cos(delta_nu));
    l = r_1 + r_2;
    m = r_1 * r_2 * (1 + cos(delta_nu));

    % Initial guesses
    p_i = k / (l + sqrt(2 * m));
    p_ii = k / (l - sqrt(2 * m));
    p_0 = (p_i + p_ii) / 2;
    
    % Calculate initial function value
    [f_p_0, ~] = tof_gauss_lambert(p_0, r_1_vec, r_2_vec, dt, mu, orbit);

    % Set up initial values based on function sign
    if f_p_0 < 0 
        p_1 = p_0;
        p_2 = p_ii;
        f_p_1 = f_p_0;
        [f_p_2, ~] = tof_gauss_lambert(p_2, r_1_vec, r_2_vec, dt, mu, orbit);
    else
        p_2 = p_0;
        p_1 = p_i;
        [f_p_1, ~] = tof_gauss_lambert(p_1, r_1_vec, r_2_vec, dt, mu, orbit);
        f_p_2 = f_p_0;
    end

    % Set tolerance, maximum iterations
    tol = 1e-12;
    iter = 1;
    n_max = 1000;
    
    % Root solving based on the selected method
    switch root_solver
        case 'linear'
            % Linear Method
            while abs(f_p_2) > tol && iter < n_max
                p_3 = p_2 - f_p_2 * (p_2 - p_1) / (f_p_2 - f_p_1);
                p_1 = p_2;
                p_2 = p_3;
                f_p_1 = f_p_2;
                [f_p_2, ~] = tof_gauss_lambert(p_2, r_1_vec, r_2_vec, dt, mu, orbit);
                iter = iter + 1;
                if iter == n_max
                    error('Maximum number of iterations reached')
                end
            end
            p = p_2;

        case 'bisection'
            % Bisection Method
            p_3 = (p_1 + p_2) / 2;
            [f_p_3, ~] = tof_gauss_lambert(p_3, r_1_vec, r_2_vec, dt, mu, orbit);
            while abs(f_p_3) > tol && iter < n_max
                if f_p_1 * f_p_3 > 0
                    p_1 = p_3;
                    f_p_1 = f_p_3;
                else
                    p_2 = p_3;
                    f_p_2 = f_p_3;
                end
                p_3 = (p_1 + p_2) / 2;
                [f_p_3, ~] = tof_gauss_lambert(p_3, r_1_vec, r_2_vec, dt, mu, orbit);
                iter = iter + 1;
                if iter == n_max
                    error('Maximum number of iterations reached')
                end
            end
            p = p_3;

        otherwise
            error('Unknown root solver method.')
    end

    % Calculate orbit parameters using solved value of p
    a = (m * k * p) / ((2 * m - l^2) * p^2 + 2 * k * l * p - k^2);
    f = 1 - (r_2 / p) * (1 - cos(delta_nu));
    g = r_1 * r_2 * sin(delta_nu) / sqrt(mu * p);
    delta_E = acos(1 - (r_1 / a) * (1 - f));
    f_dot = -sqrt(mu * a) * sin(delta_E) / (r_1 * r_2);
    g_dot = 1 - (a / r_2) * (1 - cos(delta_E));

    % Calculate initial and final velocity vectors using f and g methods
    v_1_vec = (r_2_vec - f * r_1_vec) / g;
    v_2_vec = (g_dot * r_2_vec - r_1_vec) / g;

end
