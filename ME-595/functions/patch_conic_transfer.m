function [dv_1, dv_2, dv_t] = patch_conic_transfer(mu_sun, mu_earth, mu_mars, R_earth_sun, R_mars_sun, r_p_earth, r_p_mars, a_t)
    %
    % DESCRIPTION
    %   Calculate the delta-Vs for a mission from Earth to Mars 
    %   given the semi-major axis of the transfer orbit.
    %
    % INPUTS         size    Type
    %   mu_sun      (1,1)   Double  Gravitational parameter of the Sun [DU^3/TU^2]
    %   mu_earth    (1,1)   Double  Gravitational parameter of Earth [DU^3/TU^2]
    %   mu_mars     (1,1)   Double  Gravitational parameter of Mars [DU^3/TU^2]
    %   R_earth_sun (1,1)   Double  Distance from the Sun to Earth [DU]
    %   R_mars_sun  (1,1)   Double  Distance from the Sun to Mars [DU]
    %   r_p_earth   (1,1)   Double  Radius of the parking orbit around Earth [DU]
    %   r_p_mars    (1,1)   Double  Radius of the parking orbit around Mars [DU]
    %   a_t         (1,1)   Double  Semi-major axis of the transfer orbit [DU]
    %
    % OUTPUTS            size    Type
    %   dv_1         (1,1)   Double  Delta-V required for departure from Earth [DU/TU]
    %   dv_2         (1,1)   Double  Delta-V required for capture into Mars [DU/TU]
    %   total_dv     (1,1)   Double  Total Delta-V for the mission [DU/TU]
    %
    % NOTES
    %
    % FUNCTION

    % Energy of the transfer orbit (J/kg)
    E_t = -mu_sun / (2 * a_t);  % Specific orbital energy for the transfer

    % Velocities at Earth (departure) and Mars (arrival) in the Sun-centered frame (DU/TU)
    v_1 = sqrt(2 * (mu_sun / R_earth_sun + E_t));  % Velocity at Earth's orbit in transfer
    v_earth_sun = sqrt(mu_sun / R_earth_sun);      % Circular orbit velocity of Earth around Sun
    v_inf_earth = abs(v_1 - v_earth_sun);          % Hyperbolic excess velocity at Earth (DU/TU)

    v_2 = sqrt(2 * (mu_sun / R_mars_sun + E_t));   % Velocity at Mars' orbit in transfer
    v_mars_sun = sqrt(mu_sun / R_mars_sun);        % Circular orbit velocity of Mars around Sun
    v_inf_mars = abs(v_2 - v_mars_sun);            % Hyperbolic excess velocity at Mars (DU/TU)

    % Hyperbolic escape trajectory at Earth
    E_inf_earth = v_inf_earth^2 / 2;               % Specific energy in hyperbolic escape trajectory at Earth
    v_0_earth = sqrt(2 * (mu_earth / r_p_earth + E_inf_earth));  % Velocity at perigee in the hyperbolic orbit around Earth
    h_hyp = sqrt(r_p_earth * v_inf_earth);   % Angular momentum of the hyperbolic escape
    e_hyp = sqrt(1 + (2 * E_inf_earth * h_hyp^2) / mu_earth^2);  % Eccentricity of the hyperbolic escape
    n_earth = acos(-1 / e_hyp);              % True anomaly for escape (r = infinity)
    
    % Velocity required in low Earth orbit (LEO) for escape
    v_c_earth = sqrt(mu_earth / r_p_earth);        % Circular velocity in the parking orbit (LEO)
    dv_1 = abs(v_0_earth - v_c_earth);             % Delta-V required for escape from Earth's parking orbit (DU/TU)

    % Hyperbolic capture at Mars (similar steps as Earth)
    E_inf_mars = v_inf_mars^2 / 2;                 % Specific energy in hyperbolic approach to Mars
    v_0_mars = sqrt(2 * (mu_mars / r_p_mars + E_inf_mars));  % Velocity at perigee in the hyperbolic orbit around Mars
    v_c_mars = sqrt(mu_mars / r_p_mars);           % Circular velocity in Mars' parking orbit
    dv_2 = abs(v_0_mars - v_c_mars);               % Delta-V for capture into Mars' parking orbit (DU/TU)

    % Calculate total Delta-V for the mission
    dv_t = dv_1 + dv_2;
    
end
