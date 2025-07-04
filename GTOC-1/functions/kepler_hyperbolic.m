function F = kepler_hyperbolic(M, e)
    % DESCRIPTION
    %   Solves Kepler's equation for a hyperbolic orbit using Newton's method.
    %   The function iterates to find the hyperbolic anomaly (F) that satisfies
    %   the equation e * sinh(F) - F = M, where M is the mean anomaly and e is
    %   the eccentricity of the hyperbolic orbit.
    %
    % INPUTS        size    Type    Description
    %   M           (1,1)   Double  Mean anomaly of the hyperbolic orbit
    %   e           (1,1)   Double  Eccentricity of the hyperbolic orbit
    %
    % OUTPUTS          size    Type    Description
    %   F           (1,1)   Double  Hyperbolic anomaly
    %
    % FUNCTION
    %   The function iteratively solves for the hyperbolic anomaly F using 
    %   Newton's method, starting from an initial guess based on M. The 
    %   iteration stops when the function value converges within a specified 
    %   tolerance or a maximum number of iterations is reached.

    MAXITER = 100; % maximum number of iterations
    TOL = 1e-11;   % tolerance for convergence
    Fn = (6 * M)^(1/3); % initial guess
    fn = e * sinh(Fn) - Fn - M; % function value
    iter = 1; % iteration counter
    
    while abs(fn) > TOL
        % Update guess using Newton's method
        Fn = Fn - fn / (e * cosh(Fn) - 1);
        % Update function value
        fn = e * sinh(Fn) - Fn - M;
        % Increment iteration counter
        iter = iter + 1;
        % Check for maximum iterations
        if iter >= MAXITER
            error('Did not converge');
        end
    end
    
    % Return the solution
    F = Fn;
end
