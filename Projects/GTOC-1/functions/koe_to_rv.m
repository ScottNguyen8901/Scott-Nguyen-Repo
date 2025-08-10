function [state_PQW, state_IJK] = koe_to_rv(koe, mu)
    %
    % DESCRIPTION
    %   Calculate the rotation matrix from PQW to IJK frame and the inverse matrix from IJK to PQW.
    %   Return the state (position and velocity) in both PQW and IJK frames as 6x1 matrices.
    %
    % INPUTS           size     Type
    %   a              (1,1)    Double  Semi-major axis   [DU]
    %   e              (1,1)    Double  Eccentricity      []
    %   i              (1,1)    Double  Inclination       [rad]
    %   W              (1,1)    Double  RAAN              [rad]
    %   w              (1,1)    Double  Arg of periapsis  [rad]
    %   f              (1,1)    Double  True anomaly      [rad]
    %   mu             (1,1)    Double  Gravitational parameter [DU^3/TU^2]
    %
    % OUTPUTS
    %   R_PQW_IJK      (3,3)    Double  Rotation matrix from PQW to IJK []
    %   R_IJK_PQW      (3,3)    Double  Rotation matrix from IJK to PQW []
    %   state_PQW      (6,1)    Double  State vector (position and velocity) in PQW [DU, DU/TU]
    %   state_IJK      (6,1)    Double  State vector (position and velocity) in IJK [DU, DU/TU]
    %
    % FUNCTION 

    a = koe.a;
    e = koe.e;
    i = koe.i;
    W = koe.W;
    w = koe.w;
    f = koe.f;

    R_PQW_IJK = [cos(W)*cos(w) - sin(W)*sin(w)*cos(i), -cos(W)*sin(w) - sin(W)*cos(w)*cos(i), sin(W)*sin(i);
                 sin(W)*cos(w) + cos(W)*sin(w)*cos(i), -sin(W)*sin(w) + cos(W)*cos(w)*cos(i), -cos(W)*sin(i);
                 sin(w)*sin(i),                         cos(w)*sin(i),                        cos(i)];
    
    % The rotation matrix from IJK to PQW is the transpose of R_PQW_IJK
    R_IJK_PQW = R_PQW_IJK.';

    % Handle elliptical and hyperbolic orbits
    if e < 1
        % Elliptical orbit: Use semi-major axis (a > 0)
        p = a * (1 - e^2);  % Semi-latus rectum
        r = p / (1 + e * cos(f));  % Distance from central body

        % Position in PQW frame
        r_PQW = r * [cos(f); sin(f); 0];
        
        % Velocity in PQW frame
        v_PQW = sqrt(mu / p) * [-sin(f); e + cos(f); 0];
    
    elseif e > 1
        % Hyperbolic orbit: Use semi-major axis (a < 0)
        p = -a * (e^2 - 1);  % Semi-latus rectum (note a is negative for hyperbolic orbits)
        r = p / (e * cosh(f) - 1);  % Distance from central body

        % Position in PQW frame
        r_PQW = r * [cosh(f) - e; sinh(f); 0];
        
        % Velocity in PQW frame
        v_PQW = sqrt(mu / -a) * [sinh(f); e - cosh(f); 0];
    end

    % Convert position and velocity to the IJK frame
    r_IJK = R_PQW_IJK * r_PQW;
    v_IJK = R_PQW_IJK * v_PQW;

    % Concatenate position and velocity vectors into 6x1 state vectors
    state_PQW = [r_PQW; v_PQW];
    state_IJK = [r_IJK; v_IJK];

end
