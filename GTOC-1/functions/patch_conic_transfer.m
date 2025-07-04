function [dv_p1, dv_p2, dv_t, v_inf_p1, v_inf_p2, dv_t_hel] = patch_conic_transfer(mu, mu_p1, mu_p2, R_p1, R_p2, r_p1, r_p2, a, e)
    %
    % DESCRIPTION
    %   Calculate the delta-Vs for a mission from one planet to another 
    %   given the altitudes of the parking orbits and the semi-major axes.
    %
    % INPUTS         size      Type
    %   mu          (1,1)     Double  Gravitational parameter of the Sun [DU^3/TU^2]
    %   mu_p1       (1,1)     Double  Gravitational parameter of the first planet [DU^3/TU^2]
    %   mu_p2       (1,1)     Double  Gravitational parameter of the second planet [DU^3/TU^2]
    %   R_p1        (1,1)     Double  Distance from the Sun to the first planet [DU]
    %   R_p2        (1,1)     Double  Distance from the Sun to the second planet [DU]
    %   r_p1        (1,1)     Double  Radius of the parking orbit around the first planet [DU]
    %   r_p2        (1,1)     Double  Radius of the parking orbit around the second planet [DU]
    %   a           (n,1)     Double  Vector of semi-major axes of the transfer orbits [DU]
    %   e           (n,1)     Double  Vector of eccentricities of the transfer orbits []
    %
    % OUTPUTS            size    Type
    %   dv_p1       (n,1)     Double  Delta-V required for departure from the first planet [DU/TU]
    %   dv_p2       (n,1)     Double  Delta-V required for capture into the second planet [DU/TU]
    %   dv_t        (n,1)     Double  Total Delta-V for the mission [DU/TU]
    %   v_inf_p1    (n,1)     Double  Hyperbolic excess velocity at the first planet [DU/TU]
    %   v_inf_p2    (n,1)     Double  Hyperbolic excess velocity at the second planet [DU/TU]
    %   dv_t_hel    (n,1)     Double  Total heliocentric Delta-V for the mission [DU/TU]
    %

    % Calculate the velocity at departure from the first planet
    v_1 = sqrt(2 * mu * (1 / R_p1 - 1 ./ (2 * a)));  
    % Calculate the velocity at arrival at the second planet
    v_2 = sqrt(2 * mu * (1 / R_p2 - 1 ./ (2 * a)));  
    
    % Calculate circular orbital velocities at both planets
    v_p1 = sqrt(mu / R_p1);  
    v_p2 = sqrt(mu / R_p2);  
    
    % Calculate specific angular momentum for the transfer orbits
    h = sqrt(mu .* a .* (1 - e.^2));  
    % Calculate cosine of the transfer angle using angular momentum and velocity
    cos_phi = h ./ (R_p2 .* v_2);  
    % Calculate hyperbolic excess velocity at the first planet
    v_inf_p1 = v_1 - v_p1;  
    % Calculate circular orbital velocity at the first planet for the parking orbit
    v_c_p1 = sqrt(mu_p1 / r_p1);  
    
    % Calculate specific energy at the first planet for the transfer orbit
    eps_p1 = v_inf_p1 .^ 2 / 2;  
    % Calculate velocity required for the transfer orbit at the first planet
    v_0_p1 = sqrt(2 * (mu_p1 / r_p1 + eps_p1));  
    % Calculate Delta-V required for departure from the first planet
    dv_p1 = v_0_p1 - v_c_p1;  
    
    % Calculate hyperbolic excess velocity at the second planet
    v_inf_p2 = sqrt(v_p2^2 + v_2.^2 - 2 * v_p2 .* v_2 .* cos_phi);  
    % Calculate specific energy at the second planet for the transfer orbit
    eps_p2 = v_inf_p2 .^ 2 / 2;  
    
    % Calculate velocity required for the transfer orbit at the second planet
    v_p_p2 = sqrt(2 * (mu_p2 / r_p2 + eps_p2));  
    % Calculate circular orbital velocity at the second planet for the parking orbit
    v_c_p2 = sqrt(mu_p2 / r_p2);  
    % Calculate Delta-V required for capture into the second planet
    dv_p2 = v_p_p2 - v_c_p2;  
    
    % Calculate total Delta-V for the mission
    dv_t = dv_p1 + dv_p2;  
    % Calculate total heliocentric Delta-V for the mission
    dv_t_hel = v_inf_p1 + v_inf_p2;  
end
