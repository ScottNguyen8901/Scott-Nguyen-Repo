function anomaly = kepler_solver(M, e, p, dt, mu)
    %
    % DESCRIPTION
    %   Solves Kepler's equation for elliptical, hyperbolic, and parabolic orbits 
    %   using Newton's method to find the respective orbital anomalies. This function 
    %   determines the appropriate method for solving Kepler's equation based on the 
    %   eccentricity and calculates the anomaly (eccentric, hyperbolic, or parabolic).
    %
    % INPUTS         size         Type     Description                               Units
    %   M             (1,1)        Double   Mean anomaly of the orbit                [rad]  
    %   e             (1,1)        Double   Eccentricity of the orbit                [unitless]
    %   p             (1,1)        Double   Semi-latus rectum for parabolic orbit    [DU] (optional)
    %   dt            (1,1)        Double   Time of flight for parabolic orbit       [TU] (optional)
    %   mu            (1,1)        Double   Gravitational parameter (for parabolic orbit) [DU^3/TU^2] (optional)
    %
    % OUTPUTS        size         Type     Description                               Units
    %   anomaly      (1,1)        Double   Orbital anomaly (eccentric, hyperbolic, or parabolic) [rad]
    %
    % NOTES
    %   This function calculates the appropriate anomaly based on the orbit's eccentricity:
    %   - For elliptical orbits, it solves Kepler's equation for the eccentric anomaly.
    %   - For hyperbolic orbits, it solves for the hyperbolic anomaly.
    %   - For parabolic orbits, it solves for the parabolic anomaly using the given 
    %     semi-latus rectum, time of flight, and gravitational parameter.
    %   The function uses Newton's method to iteratively solve the respective equation 
    %   and terminates once the error is below a specified tolerance.
    %

    if nargin < 4, p = 0; dt = 0; mu = 0; end
    
    tol = 1e-11; % Common tolerance for all cases
    max_iter = 100; % Common max iterations for all cases
    
    if e < 1  % Elliptical Orbit
        E_n = M; % Initial guess
        for iter = 1:max_iter
            f_n = E_n - e*sin(E_n) - M;
            if abs(f_n) < tol, anomaly = E_n; return; end
            E_n = E_n - f_n / (1 - e * cos(E_n));
        end
        error('Elliptical orbit did not converge');
        
    elseif e > 1  % Hyperbolic Orbit
        Fn = (6 * M)^(1/3); % Initial guess
        for iter = 1:max_iter
            fn = e * sinh(Fn) - Fn - M;
            if abs(fn) < tol, anomaly = Fn; return; end
            Fn = Fn - fn / (e * cosh(Fn) - 1);
        end
        error('Hyperbolic orbit did not converge');
        
    else  % Parabolic Orbit (e = 1)
        if p == 0 || dt == 0 || mu == 0, error('For parabolic orbit, p, dt, and mu must be provided'); end
        Dn = M; % Initial guess
        for iter = 1:max_iter
            fn = (p * Dn + Dn^3 / 3) - 2 * sqrt(mu) * dt;
            if abs(fn) < tol, anomaly = Dn; return; end
            Dn = Dn - fn / (p + Dn^2);
        end
        error('Parabolic orbit did not converge');
    end
end
