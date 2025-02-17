function state_vectors = orbit_prop_fg(r_1_vec, v_1_vec, time_vec, mu)
    %
    % DESCRIPTION
    %   Propagate an orbit using position and velocity vectors at an initial point
    %   and a time vector. The function computes the position and velocity vectors
    %   at each time in the time vector using the orbital elements and the f and g method.
    %
    % INPUTS            size         Type    Description                   Units
    %   r_1_vec         (3,1)        Double  Initial position vector       [DU]
    %   v_1_vec         (3,1)        Double  Initial velocity vector       [DU/TU]
    %   mu              (1,1)        Double  Gravitational parameter       [DU^3/TU^2]
    %   time_vec        (N,1)        Double  Vector of times of flight     [TU]
    %
    % OUTPUTS           size         Type    Description                   Units
    %   state_vectors   (N,7)        Double  Time, Position and Velocity   [TU, DU, DU/TU]
    %
    % NOTES
    %   This function uses the Gauss f and g method for orbit propagation. It requires
    %   initial position and velocity vectors and propagates them forward in time for each
    %   time value in the time vector. Orbital elements are computed and used to obtain
    %   the position and velocity vectors for each time.
    %
    % FUNCTION

    % Step 1: Compute initial quantities
    r = norm(r_1_vec);
    v = norm(v_1_vec);
    
    % Convert position and velocity vectors to orbital elements
    [koe] = rv2koe(r_1_vec, v_1_vec, mu);
    
    % #3: Semimajor axis
    a = koe.a;
    n = sqrt(mu / a^3);  % Mean motion
    
    % #4: Eccentricity
    e = koe.e;
    p = a * (1 - e^2);  % Semi-latus rectum
    
    % #5: True anomaly at the initial state
    nu = koe.nu;
    
    % #6: Time since periapsis passage (not an output in this case)
    cosE = (a * e + r * cos(nu)) / a;
    E = acos(cosE);  % Eccentric anomaly
    M = E - e * sin(E);  % Mean anomaly
    t_1 = M / n ; % Time since periapsis
    
    % Pre-allocate the state vector array
    N = length(time_vec);
    state_vectors = zeros(N, 7);  % [time, r, v]
    
    % Loop over the time vector
    for i = 1:N
        % Calculate the total time after periapsis passage
        t_2 = t_1 + time_vec(i);  % Total time after periapsis passage
        M_2 = n * t_2;  % New mean anomaly after the time of flight
        E_2 = kepler(M_2, e);  % Solve Kepler's equation for new eccentric anomaly
        r_2 = a * (1 - e * cos(E_2));  % New radius at time of flight
        cos_nu_2 = (e - cos(E_2)) / (e * cos(E_2) - 1);  % Cosine of new true anomaly
        sin_nu_2 = a * sqrt(1 - e^2) * sin(E_2) / r_2;  % Sine of new true anomaly
        nu_2 = atan2(sin_nu_2, cos_nu_2);  % Compute new true anomaly
        
        % Ensure the true anomaly is positive
        if nu_2 < 0
            nu_2 = 2 * pi + nu_2;
        end
        
        % #8: Value of f after the time of flight
        dnu = nu_2 - nu;  % Change in true anomaly
        f = 1 - r_2 / p * (1 - cos(dnu));  % f function
        
        % #9: New position and velocity vectors after the time of flight
        g = r_2 * r * sin(dnu) / sqrt(mu * p);  % g function
        fdot = sqrt(mu / p) * tan(dnu / 2) * ((1 - cos(dnu)) / p - 1 / r - 1 / r_2);  % f-dot
        gdot = 1 - r / p * (1 - cos(dnu));  % g-dot
        
        % New position and velocity vectors
        r_2_vec = f * r_1_vec + g * v_1_vec;
        v_2_vec = fdot * r_1_vec + gdot * v_1_vec;
        
        % Store the results (time, position, and velocity)
        state_vectors(i, :) = [time_vec(i), r_2_vec', v_2_vec'];
    end
end
