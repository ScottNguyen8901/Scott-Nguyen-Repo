function dydt = two_body_ode_drag(t, y, rho_table)
    % DESCRIPTION
    %   Computes the derivative of the state vector for two-body motion 
    %   including atmospheric drag based on an exponential density model.
    %
    % INPUTS       size     Type       Description                              Units
    %   t          (1,1)    Double     Time                                     [s]
    %   y          (6,1)    Double     State vector: [position; velocity]        [km; km/s]
    %   rho_table  (M,2)    Double     Precomputed density values at altitudes  [kg/km^3]
    %
    % OUTPUTS      size     Type       Description                              Units
    %   dydt       (6,1)    Double     Derivative of the state vector:          [km/s; km/s^2]
    %                                  [velocity; acceleration]                   
    %
    % FUNCTION
    %   The function computes the time derivatives of the satellite's position 
    %   and velocity considering two-body gravitational motion and atmospheric drag.
    %   The atmospheric drag is modeled using a precomputed density table, which 
    %   significantly improves computational efficiency.

    constants;

    % Extract position and velocity vectors
    r_vec = y(1:3);           % Position vector [km]
    v_vec = y(4:6);           % Velocity vector [km/s]
    r = norm(r_vec);          % Magnitude of position vector [km]

    % Gravitational acceleration
    a_gravity = -mu_E * r_vec / r^3;

    % Compute altitude above Earth's surface [km]
    h = r - R_E;

    % Find the index corresponding to the closest altitude in the rho_table
    % Since the altitude is represented by the row index in rho_table,
    % we assume that the altitude (h) matches the row number.
    altitude_index = round(h);  % Round the altitude to the nearest integer (index)
    
    % Ensure the altitude_index does not exceed the bounds of the rho_table
    if altitude_index > size(rho_table, 1)
        altitude_index = size(rho_table, 1);  % Set to maximum if out of bounds
    elseif altitude_index < 1
        altitude_index = 1;  % Set to minimum if out of bounds
    end
    
    % Get the atmospheric density for the corresponding altitude
    rho = rho_table(altitude_index, 2);  % Get density at closest altitude [kg/km^3]

    % Relative velocity (satellite velocity - Earth's rotation x position)
    v_rel = v_vec - cross(omega_E_vec, r_vec);  % [km/s]
    v_rel_mag = norm(v_rel);

    % Drag acceleration
    if rho == 0
        a_drag = [0; 0; 0];  % No atmospheric drag if density is zero
    else
        a_drag = -0.5 * Cd * A / m * rho * v_rel_mag * v_rel * 1e3;  % [km/s^2] Drag force
    end

    % Total acceleration
    a_total = a_gravity + a_drag;

    % Derivative of state
    dydt = [v_vec; a_total];  % Return velocity and acceleration
end