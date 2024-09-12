function [dv1, dv2, dvt, T_transfer] = hohmann_transfer(r1, r2, mu)
    %
    % DESCRIPTION
    %   Calculate the delta-Vs and time of flight for a Hohmann transfer between 
    %   two circular orbits.
    %
    % INPUTS    size    Type
    %   r1      (1,1)   Double  Radius of initial orbit [DU]
    %   r2      (1,1)   Double  Radius of final orbit   [DU]
    %   mu      (1,1)   Double  Gravitational parameter [DU^3/TU^2]
    %
    % OUTPUTS       size    Type
    %   dv1         (1,1)   Double  First delta-V  [DU/TU]
    %   dv2         (1,1)   Double  Second delta-V [DU/TU]
    %   dvt         (1,1)   Double  Total delta-V  [DU/TU]
    %   T_transfer  (1,1)   Double  Time of flight [TU]
    %

    % Velocities in initial and final orbits
    v1 = sqrt(mu / r1); % velocity in the initial orbit
    v2 = sqrt(mu / r2); % velocity in the final orbit

    % Semi-major axis of the transfer orbit
    a_transfer = (r1 + r2) / 2;

    % Velocities at periapsis and apoapsis of the transfer orbit
    v_periapsis = sqrt(mu * (2 / r1 - 1 / a_transfer)); % velocity at periapsis (r1
    v_apoapsis = sqrt(mu * (2 / r2 - 1 / a_transfer)); % velocity at apoapsis (r2)

    % Delta-V calculations
    dv1 = v_periapsis - v1; % first burn (at r1)
    dv2 = v2 - v_apoapsis;  % second burn (at r2)
    dvt = abs(dv1) + abs(dv2); % total Delta-V for the Hohmann transfer

    % Time of Flight (half the orbital period of the transfer orbit)
    T_transfer = pi * sqrt(a_transfer^3 / mu); % time for half the transfer orbit [s]

end
