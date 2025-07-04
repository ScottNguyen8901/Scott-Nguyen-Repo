classdef orb_ele
    % Class to represent a celestial body's orbital elements

    properties
        a   % Semi-major axis (a) in Astronomical Units (AU)
        e   % Orbital eccentricity (e)
        i   % Inclination (i) in degrees
        w  % Argument of pericenter (ω) in degrees
        W  % Right Ascension of the Ascending Node (Ω) in degrees
        M     % Mean anomaly (M) at epoch in degrees
        JD    % Initial Julian Date
        f     % True anomaly (f) in degrees
    end

    methods
        % Constructor to initialize the orbital elements
        function obj = orb_ele(a, e, i, omega, RAAN, M, JD)
            if nargin > 0
                obj.a = a;
                obj.e = e;
                obj.i = i;
                obj.w = omega;
                obj.W = RAAN;
                obj.M = M;
                obj.JD = JD;
                % Calculate the true anomaly using the provided code
                F = kepler_hyperbolic(obj.M, obj.e);
                obj.f = ecc_to_true_anom(F, obj.e);  % True anomaly calculation
            end
        end

        % Method to display orbital elements
        function displayInfo(obj)
            fprintf('Orbital Elements:\n');
            fprintf('  Semi-major axis (sma): %.6f km\n', obj.a);
            fprintf('  Eccentricity (ecc): %.6f\n', obj.e);
            fprintf('  Inclination (inc): %.5f rad.\n', obj.i);
            fprintf('  Argument of pericenter (argp): %.5f rad.\n', obj.w);
            fprintf('  RAAN (raan): %.5f rad.\n', obj.W);
            fprintf('  Mean anomaly (M): %.5f rad.\n', obj.M);
            fprintf('  True anomaly (f): %.5f rad.\n', obj.f);  % Display true anomaly
            fprintf('  Initial Julian Date (JD): %.0f\n', obj.JD);
        end
    end
end
