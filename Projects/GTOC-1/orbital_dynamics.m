function dydt = orbital_dynamics(t, state, mu)
    % Extract position and velocity from state vector
    r = state(1:3);  % Position
    v = state(4:6);  % Velocity
    
    % Compute the gravitational force (acceleration)
    r_norm = norm(r);  % Distance from the Sun
    a = -mu * r / r_norm^3;  % Gravitational acceleration
    
    % Return the derivatives: velocity and acceleration
    dydt = [v; a];  % [dr/dt; dv/dt]
end
