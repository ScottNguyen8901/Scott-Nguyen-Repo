function [dv_1, dv_2, dv_t] = patch_conic_transfer(mu_sun, mu_earth, mu_mars, R_sun_earth, R_sun_mars, r_p_earth, r_p_mars, p, e)
    %
    % DESCRIPTION
    %   Calculate the delta-Vs for a mission from Earth to Mars 
    %   given the semi-major axis of the transfer orbit.
    %
    % INPUTS         size    Type
    %   mu_sun      (1,1)   Double  Gravitational parameter of the Sun [DU^3/TU^2]
    %   mu_earth    (1,1)   Double  Gravitational parameter of Earth   [DU^3/TU^2]
    %   mu_mars     (1,1)   Double  Gravitational parameter of Mars    [DU^3/TU^2]
    %   R_earth_sun (1,1)   Double  Distance from the Sun to Earth     [DU]
    %   R_mars_sun  (1,1)   Double  Distance from the Sun to Mars      [DU]
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
    a_t = p / (1 - e^2);        % Semi-major axis of transfer orbit
    E_t = -mu_sun / (2 * a_t);  % Specific orbital energy for the transfer

    % Velocities at Earth (departure) and Mars (arrival) in the Sun-centered frame (DU/TU)
    v_1 = sqrt(2 * (mu_sun / R_sun_earth + E_t));  % Velocity at Earth's orbit in transfer
    v_sun_earth = sqrt(mu_sun / R_sun_earth);      % Circular orbit velocity of Earth around Sun
    v_inf_earth = abs(v_1 - v_sun_earth);          % Hyperbolic excess velocity at Earth (DU/TU)

    v_2 = sqrt(2 * (mu_sun / R_sun_mars + E_t));   % Velocity at Mars' orbit in transfer
    v_sun_mars = sqrt(mu_sun / R_sun_mars);        % Circular orbit velocity of Mars around Sun
    v_inf_mars = abs(v_2 - v_sun_mars);            % Hyperbolic excess velocity at Mars (DU/TU)
    
    % Calculate sphere of influence
    R_SOI_earth = R_sun_earth * (mu_earth / mu_sun) ^ (2/5);
    R_SOI_mars  = R_sun_mars * (mu_mars / mu_sun) ^ (2/5);

    % Hyperbolic escape trajectory at Earth
    E_inf_earth = v_inf_earth^2 / 2 - mu_earth / R_SOI_earth;               % Specific energy in hyperbolic escape trajectory at Earth
    v_0_earth = sqrt(2 * (mu_earth / r_p_earth + E_inf_earth));  % Velocity at perigee in the hyperbolic orbit around Earth
    h_hyp_earth = r_p_earth * v_0_earth;   % Angular momentum of the hyperbolic escape
    e_hyp_earth = sqrt(1 + (2 * E_inf_earth * h_hyp_earth^2) / mu_earth^2);  % Eccentricity of the hyperbolic escape
    nu_earth = acos(-1 / e_hyp_earth);              % True anomaly for escape (r = infinity)
    
    % Velocity required in low Earth orbit (LEO) for escape
    v_c_earth = sqrt(mu_earth / r_p_earth);        % Circular velocity in the parking orbit (LEO)
    dv_1 = abs(v_0_earth - v_c_earth);             % Delta-V required for escape from Earth's parking orbit (DU/TU)

    % Hyperbolic capture at Mars (similar steps as Earth)
    E_inf_mars = v_inf_mars^2 / 2 - mu_mars / R_SOI_mars;                 % Specific energy in hyperbolic approach to Mars
    v_0_mars = sqrt(2 * (mu_mars / r_p_mars + E_inf_mars));  % Velocity at perigee in the hyperbolic orbit around Mars
    h_hyp_mars = r_p_mars * v_0_mars;   % Angular momentum of the hyperbolic escape
    e_hyp_mars = sqrt(1 + (2 * E_inf_mars * h_hyp_mars^2) / mu_mars^2);  % Eccentricity of the hyperbolic escape
    nu_mars = acos(-1 / e_hyp_mars);              % True anomaly for escape (r = infinity)

    v_c_mars = sqrt(mu_mars / r_p_mars);           % Circular velocity in Mars' parking orbit
    dv_2 = abs(v_0_mars - v_c_mars);               % Delta-V for capture into Mars' parking orbit (DU/TU)

    % Calculate total Delta-V for the mission
    dv_t = dv_1 + dv_2;
    
end