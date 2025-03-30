function dydt = two_body_ode(t, y, mu)
    % TWO_BODY_ODE Computes the derivative for two-body motion
    % Inputs:
    %   y    - State vector [position; velocity]
    %   mu   - Gravitational parameter (GM)
    % Output:
    %   dydt - Derivative of the state vector

    r = y(1:3);  % Position vector
    v = y(4:6);  % Velocity vector
    
    % Compute the acceleration due to gravity (Newton's law of gravitation)
    r_norm = norm(r);  % Magnitude of position vector
    a = -mu * r / r_norm^3;  % Acceleration vector (Newton's law)
    
    % Derivative of position is velocity
    % Derivative of velocity is acceleration
    dydt = [v; a];
end