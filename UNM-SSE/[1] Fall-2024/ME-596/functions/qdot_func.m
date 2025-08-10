function q_dot = qdot_func(t, q, w)
    % DESCRIPTION
    %   Computes the time derivative of a quaternion for integration using
    %   the equation of motion.
    %
    % INPUTS       size    Type    Description
    %   t          (1, 1)   Double  Current time (required for ode45 but not used in this function)
    %   q          (4, 1)   Double  Quaternion in the form [qx, qy, qz, qw]'
    %   w          (3, 1)   Double  Angular velocity vector in the body frame
    %
    % OUTPUTS      size    Type    Description
    %   qdot       (4, 1)   Double  Time derivative of the quaternion [dq_x; dq_y; dq_z; dq_w]'
    %
    % FUNCTION
    
    % Quaternion derivative calculation
    q_dot = [skew(q(1:3)) + q(4) * eye(3); -q(1:3)'] * w;
end
