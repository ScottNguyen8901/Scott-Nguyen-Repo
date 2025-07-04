function [r, v] = koe_to_rv(koe, mu)
    %
    % DESCRIPTION
    %   Converts orbital elements (semi-major axis, eccentricity, inclination, right 
    %   ascension of ascending node, argument of periapsis, and true anomaly) into 
    %   the position and velocity vectors (state vector). 
    %   This function handles elliptical, parabolic, and hyperbolic orbits.
    %
    % INPUTS         size         Type     Description                               Units
    %   koe          (1,1)        Struct   Structure containing the Keplerian orbital
    %                                     elements with the following fields:
    %                                     - a: Semi-major axis                       [DU]    
    %                                     - e: Eccentricity                         
    %                                     - i: Inclination
    %                                     - w: Argument of Periapsis
    %                                     - W: Right Ascension of Ascending Node
    %                                     - f: True Anomaly
    %   mu           (1,1)        Double   Gravitational Parameter                   [DU^3/TU^2]
    %
    % OUTPUTS        size         Type     Description                               Units
    %   r            (3,1)        Double   Position vector                           [DU]
    %   v            (3,1)        Double   Velocity vector                           [DU/TU]
    %
    % NOTES
    %   This function calculates the position and velocity vectors for a given orbit 
    %   based on the orbital elements. The conversion is done for elliptical, parabolic, 
    %   and hyperbolic orbits. It uses the appropriate formulae for each case and then 
    %   transforms the results from the perifocal coordinate system to the ECI frame 
    %   using a rotation matrix.
    %

    % Extract orbital elements from the koe structure
    a = koe.a;
    e = koe.e;
    i = koe.i;
    w = koe.w;
    W = koe.W;
    f = koe.f;
    
    % Convert angles from degrees to radians
    i = deg2rad(i);
    W = deg2rad(W);
    w = deg2rad(w);
    f = deg2rad(f);
    
    % Compute the orbital radius (r) and velocity magnitude (v)
    if e < 1
        % Elliptical orbit
        p = a * (1 - e^2);  % Semi-latus rectum
        r = p / (1 + e * cos(f));  % Orbital radius
        
        % Velocity magnitude for elliptical orbits
        v = sqrt(mu * (2 / r - 1 / a));
    elseif e == 1
        % Parabolic orbit
        r = a * (1 + e * cos(f));  % For parabolic orbits, a = infinity
        
        % Parabolic velocity magnitude
        v = sqrt(mu * (2 / r));
    else
        % Hyperbolic orbit
        r = a * (e - 1) / (1 + e * cos(f));  % For hyperbolic orbits
        
        % Hyperbolic velocity magnitude
        v = sqrt(mu * (2 / r + 1 / a));
    end
    
    % True anomaly position components (in the perifocal frame)
    r_pqw = [r * cos(f); r * sin(f); 0];
    v_pqw = [v * (-sin(f)); v * (e + cos(f)); 0];
    
    % Rotation matrix from perifocal to ECI frame
    R = rot_mat(W, w, i);
    
    % Convert position and velocity to the ECI frame
    r_eci = R * r_pqw;
    v_eci = R * v_pqw;
    
    % Outputs
    r = r_eci;
    v = v_eci;
end