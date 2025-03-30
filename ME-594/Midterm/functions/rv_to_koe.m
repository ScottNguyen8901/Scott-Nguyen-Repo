function [koe] = rv_to_koe(r_vec, v_vec, mu)
    %
    % DESCRIPTION
    %   Convert state vector (position and velocity) into Keplerian orbital elements
    %
    % INPUTS    
    %   r_vec (3,1)   Position vector [DU]   
    %   v_vec (3,1)   Velocity vector [DU^3/TU^2]
    %   mu    (1,1)   Gravitational parameter [DU^3/TU^2]
    %
    % OUTPUTS
    %   koe     {struct}   Keplerian orbital elements structure
    %    .a     (1,1)      Semi-major axis [DU]
    %    .e     (1,1)      Eccentricity []
    %    .i     (1,1)      Inclination [rad]
    %    .W     (1,1)      RAAN [rad]
    %    .w     (1,1)      Argument of periapsis [rad]
    %    .f     (1,1)      True anomaly [rad]
    %    .fpa   (1,1)      Flight path angle [rad]
    %
    % FUNCTION IMPLEMENTATION
    
    % Unit vectors for i, j, k axes
    I = [1 0 0];
    J = [0 1 0];
    K = [0 0 1];
    
    % Check if the input vectors are non-zero
    if norm(r_vec) == 0 || norm(v_vec) == 0
        error('Position and velocity vectors must be non-zero.');
    end
    
    % Norm of r and v
    r = norm(r_vec);
    v = norm(v_vec);
    
    % Calculating semi-major axis using the Vis-Viva equation
    a = 1 / ((2 / r) - ((v^2) / mu));
    
    % Calculating eccentricity vector and its magnitude
    e_vec = (((v^2) / mu) - (1 / r)) * r_vec - (1 / mu) * dot(r_vec, v_vec) * v_vec;
    e = norm(e_vec);
    
    % Angular momentum vector
    h_vec = cross(r_vec, v_vec);
    
    % Inclination (i) – angle between h_vec and K (z-axis)
    i = acos(dot(h_vec / norm(h_vec), K));
    
    % Node vector (n) – cross product of K and h_vec
    n = cross(K, h_vec);
    
    % Right Ascension of Ascending Node (RAAN, W)
    if dot(n, J) < 0
        W = acos(dot(n / norm(n), I));
        W = 2 * pi - W;  % Adjust RAAN to correct quadrant
    else
        W = acos(dot(n / norm(n), I));
    end
    
    % Argument of Periapsis (w)
    if dot(e_vec, K) < 0
        w = acos(dot(n / norm(n), e_vec / norm(e_vec)));
        w = 2 * pi - w;  % Adjust argument of periapsis to correct quadrant
    else
        w = acos(dot(n / norm(n), e_vec / norm(e_vec)));
    end
    
    % True Anomaly (f)
    if dot(r_vec, v_vec) < 0
        f = acos(dot(r_vec / norm(r_vec), e_vec / norm(e_vec)));
        f = 2 * pi - f;  % Adjust true anomaly to correct quadrant
    else
        f = acos(dot(r_vec / norm(r_vec), e_vec / norm(e_vec)));
    end
    
    % Flight Path Angle (fpa) – angle between velocity vector and angular momentum vector
    fpa = acos(norm(h_vec) / (norm(r_vec) * norm(v_vec)));
    
    % Populating KOE output structure
    koe.a = a;     % Semi-major axis
    koe.e = e;     % Eccentricity
    koe.i = i;     % Inclination
    koe.W = W;     % RAAN
    koe.w = w;     % Argument of periapsis
    koe.f = f;     % True anomaly
    koe.fpa = fpa; % Flight path angle
end