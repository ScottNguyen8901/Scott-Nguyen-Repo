function oe = rv_to_koe(x, mu)
%
% DESCRIPTION
%   Convert state vector (position and velocity) into Keplerian orbital elements
%
% INPUTS       Size    Type      Description                    Units
%   x          (n,6)   (Double)   State vector (position and velocity)  [km, km/s]
%   mu         (1,1)   (Double)   Gravitational parameter (central body)  [km^3/s^2]
%
% OUTPUTS      Size    Type      Description                          Units
%   orbital_elements (n,6)   (Double)   Matrix of orbital elements         []
%      .a      (1,1)   (Double)   Semi-major axis                     [km]
%      .e      (1,1)   (Double)   Eccentricity                        []
%      .i      (1,1)   (Double)   Inclination                         [rad]
%      .Omega  (1,1)   (Double)   RAAN (Right Ascension of Ascending Node)  [rad]
%      .omega  (1,1)   (Double)   Argument of periapsis               [rad]
%      .nu     (1,1)   (Double)   True anomaly                        [rad]
%
% FUNCTION IMPLEMENTATION
%

    % Initialize matrix to store orbital elements for each state vector
    num_steps = size(x, 1);  % Number of time steps
    oe = zeros(num_steps, 6);  % To store [a, e, i, Omega, omega, nu]
    
    for idx = 1:num_steps
        r = x(idx, 1:3);  % Position vector (km)
        v = x(idx, 4:6);  % Velocity vector (km/s)
        
        % Magnitude of position and velocity
        r_mag = norm(r);  % Distance to the central body (km)
        v_mag = norm(v);  % Speed (km/s)
        
        % Specific angular momentum
        h = cross(r, v);  % magnitude of h is specific angular momentum (km^2/s)
        h_mag = norm(h);
        
        % Eccentricity vector
        e_vec = (cross(v, h) / mu) - (r / r_mag);
        e = norm(e_vec);  % Eccentricity (dimensionless)
        
        % Semi-major axis (a) using vis-viva equation
        a = 1 / ((2 / r_mag) - (v_mag^2 / mu));  % Semi-major axis (km)
        
        % Inclination (i)
        i = acos(h(3) / h_mag);  % Inclination (radians)
        
        % Right Ascension of Ascending Node (Omega)
        Omega = atan2(h(1), -h(2));  % RAAN (radians)
        
        % Argument of Periapsis (omega)
        N = cross([0 0 1], h);  % Node vector (km)
        if norm(N) == 0
            omega = 0;  % Periapsis argument undefined in circular orbit
        else
            omega = acos(dot(N, e_vec) / (norm(N) * e));  % Argument of periapsis (radians)
            if e_vec(3) < 0
                omega = 2 * pi - omega;  % Correct for angle direction
            end
        end
        
        % True Anomaly (nu)
        nu = acos(dot(e_vec, r) / (e * r_mag));  % True anomaly (radians)
        if dot(r, v) < 0
            nu = 2 * pi - nu;  % If the object is moving away from periapsis
        end
        
        % Store the orbital elements for this time step
        oe(idx, :) = [a, e, i, Omega, omega, nu];
    end
end