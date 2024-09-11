function omega = calcArgPer(n, e)
    %
    % DESCRIPTION
    %   Calculate the argument of perigee for an orbit given the eccentricity vector
    %
    % INPUTS    
    %   e      (3,1)   Eccentricity vector []
    %   n      (3,1)   Node line vector    []
    %
    % OUTPUTS
    %   omega  (1,1)   Argument of perigee [rad]
    %
    % NOTES
    %   The node line vector is assumed to be along the x-axis for this calculation
    %
    % FUNCTION 

    omega = acos(dot(e, n) / (norm(e) * norm(n)));
    
    % Adjust omega to the correct quadrant
    if e(3) < 0
        omega = 2 * pi - omega;  % Ensure omega is in the range [0, 2*pi]
    end
    
end