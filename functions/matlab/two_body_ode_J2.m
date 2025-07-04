function dydt = two_body_ode_J2(t, y)
%
% DESCRIPTION
%   Computes the derivative of the state vector for two-body motion
%   including the J2 perturbation due to Earth's oblateness.
%
% INPUTS   Size     Type      Description                       Units
%   t      (1,1)    (Double)   Time                             [TU]
%   y      (6,1)    (Double)   State vector: [pos; vel]         [DU, DU/TU]
%
% OUTPUTS  Size     Type      Description                       Units
%   dydt   (6,1)    (Double)   Derivative of the state vector   [DU/TU, DU/TU^2]
%
% FUNCTION IMPLEMENTATION

    constants;

    r_vec = y(1:3);           % Position vector [x; y; z]
    v_vec = y(4:6);           % Velocity vector [vx; vy; vz]
    r = norm(r_vec);          % Magnitude of position vector

    x = r_vec(1);
    y_ = r_vec(2);
    z = r_vec(3);

    % J2 perturbation acceleration vector
    factor = (3/2) * J_2 * mu_E * R_E^2 / r^5;
    zx_ratio_sq = (z/r)^2;

    a_J2_x = factor * x * (5*zx_ratio_sq - 1);
    a_J2_y = factor * y_ * (5*zx_ratio_sq - 1);
    a_J2_z = factor * z * (5*zx_ratio_sq - 3);
    a_J2 = [a_J2_x; a_J2_y; a_J2_z];

    % Two-body gravitational acceleration
    a_gravity = -mu_E * r_vec / r^3;

    % Total acceleration
    a_total = a_gravity + a_J2;

    dydt = [v_vec; a_total];

end