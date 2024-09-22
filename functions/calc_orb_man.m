function [dv_t, dt, orbit_type] = calc_orb_man(r_1, r_2, e_vec, mu)
    % Preallocate arrays for results
    dv_t = zeros(size(e_vec));
    dt = zeros(size(e_vec));
    orbit_type = zeros(size(e_vec)); % 1 for elliptical, 2 for parabolic, 3 for hyperbolic
        
    r_p = r_1;
    % Velocities in initial and final orbits
    v_1 = sqrt(mu / r_1); % velocity in the parking orbit
    v_2 = sqrt(mu / r_2); % velocity in the mission orbit

    for i = 1:length(e_vec)
        e = e_vec(i);
        a = r_p / (1 - e);
        
        if e < 1
            % Elliptical Orbits
            v1_ellip = sqrt(mu * (2 / r_1 - 1 / a)); % velocity at periapsis (r1)
            v2_ellip = sqrt(mu * (2 / r_2 - 1 / a)); % velocity at apoapsis (r2)
            dv1_ellip = v1_ellip - sqrt(mu / r_1); % first burn (at periapsis)
            dv2_ellip = sqrt(mu / r_2) - v2_ellip; % second burn (at apoapsis)
            dv_t(i) = abs(dv1_ellip) + abs(dv2_ellip);

            E = acos((1 - r_2 / a) / e); % Eccentric anomaly
            dt(i) = sqrt(a^3 / mu) * (E - e * sin(E));
            orbit_type(i) = 1; % 1 for elliptical

        elseif e == 1
            % Parabolic Orbits
            v1_par = sqrt(mu * (2 / r_1));
            v2_par = sqrt(mu * (2 / r_2));
            dv1_par = v1_par - sqrt(mu / r_1); % first burn (at periapsis)
            dv2_par = sqrt(mu / r_2) - v2_par; % second burn (at apoapsis)
            dv_t(i) = abs(dv1_par) + abs(dv2_par);

            % Semi-latus rectum for a parabolic orbit
            p = 2 * r_1; % Semi-latus rectum
            nu2_rad = acos((p - r_2) / r_2); % Calculate from orbital geometry
            D = sqrt(p) * tan(nu2_rad / 2);
            dt(i) = 1 / (2 * sqrt(mu)) * (p * D + D^3 / 3);
            orbit_type(i) = 2; % 2 for parabolic

        else
            % Velocities at periapsis and apoapsis of the elliptical transfer orbit
            v_1_hyp = sqrt(mu * (2 / r_1 - 1 / a)); % velocity at periapsis (r1)
            v_2_hyp = sqrt(mu * (2 / r_2 - 1 / a)); % velocity at apoapsis (r2)

            dv_1 = v_1_hyp - v_1; % first burn (at periapsis)
            dv_2 = v_2 - v_2_hyp;  % second burn (at apoapsis)
            dv_t(i) = abs(dv_1) + abs(dv_2); % total Delta-V for the transfer
            
            F = acosh((1 - r_2 / a) / e);
            dt(i) = sqrt(-a^3 / mu)*(e * sinh(F) - F);

            orbit_type(i) = 3; % 3 for hyperbolic
        end
    end
end
