function [r_vec, v_vec] = koe_to_rv(koe, mu)
%
% DESCRIPTION
%   Convert classical orbital elements into state vectors (position and velocity)
%   using normalized units (DU/TU).
%
% INPUTS    
%   koe   (struct) Struct containing the orbital elements:
%                 - koe.a  : Semi-major axis [DU]
%                 - koe.e  : Eccentricity []
%                 - koe.i  : Inclination [rad]
%                 - koe.W  : Right Ascension of Ascending Node (RAAN) [rad]
%                 - koe.w  : Argument of Periapsis [rad]
%                 - koe.f  : True anomaly [rad]
%   mu    (1,1)    Gravitational parameter [DU³/TU²]
%
% OUTPUTS
%   r_vec (3,1)    Position vector in ECI frame [DU]
%   v_vec (3,1)    Velocity vector in ECI frame [DU/TU]
%
% NOTES
%   All angles must be in radians.
%   Output vectors are expressed in the Earth-Centered Inertial (ECI) frame.
%

% Extract orbital elements from the struct
a = koe.a;
e = koe.e;
i = koe.i;
W = koe.W;
w = koe.w;
f = koe.f;

% Compute specific angular momentum
h = sqrt(mu * a * (1 - e^2)); % [DU²/TU]

% Position in perifocal frame [DU]
r_perifocal = (h^2 / mu) * (1 / (1 + e * cos(f))) * [cos(f); sin(f); 0];

% Velocity in perifocal frame [DU/TU]
v_perifocal = (mu / h) * [-sin(f); (e + cos(f)); 0];

% Rotation matrix from perifocal to ECI frame
R = [cos(W)*cos(w) - sin(W)*sin(w)*cos(i), ...
     -cos(W)*sin(w) - sin(W)*cos(w)*cos(i), ...
      sin(W)*sin(i);
     
     sin(W)*cos(w) + cos(W)*sin(w)*cos(i), ...
     -sin(W)*sin(w) + cos(W)*cos(w)*cos(i), ...
     -cos(W)*sin(i);
     
     sin(w)*sin(i), ...
      cos(w)*sin(i), ...
      cos(i)];

% Transform to ECI frame
r_vec = R * r_perifocal; % [DU]
v_vec = R * v_perifocal; % [DU/TU]

end