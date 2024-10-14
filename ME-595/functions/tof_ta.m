function tof = tof_ta(a, e, nu_1, nu_2, mu)
    %
    % DESCRIPTION
    %   Compute the time of flight between two true anomalies for an elliptical orbit.
    %   The function uses the eccentricity, semi-major axis, gravitational parameter, and
    %   the true anomalies (in degrees) to compute the time taken to travel between them.
    %
    % INPUTS    size     Type    Description                   Units   
    %   e       (1,1)    Double  Eccentricity                  []
    %   a       (1,1)    Double  Semi-major axis               [DU]
    %   mu      (1,1)    Double  Gravitational parameter       [DU^3/TU^2]
    %   nu_1    (1,1)    Double  True anomaly 1                [rad]
    %   nu_2    (1,1)    Double  True anomaly 2                [rad]
    %
    % OUTPUTS   size     Type    Description                   Units
    %   tof     (1,1)    Double  TOF between the two TA        [TU]
    %
    % NOTES
    %
    % FUNCTION
    
    if e > 0 && e < 1
        % Calculate the eccentric anomaly E for both true anomalies
        tan_nu_1_over_2 = tan(nu_1 / 2);
        sqrt_term = sqrt((1 + e) / (1 - e));
        tan_E1_over_2 = tan_nu_1_over_2 / sqrt_term;
        E_1 = 2 * atan(tan_E1_over_2);
    
        tan_nu_2_over_2 = tan(nu_2 / 2);
        tan_E_2_over_2 = tan_nu_2_over_2 / sqrt_term;
        E_2 = 2 * atan(tan_E_2_over_2);
    
        % Calculate the mean anomalies M from E
        M_1 = E_1 - e * sin(E_1);
        M_2 = E_2 - e * sin(E_2);
    
        % Normalize the mean anomalies M within the range [0, 2*pi]
        M_1 = mod(M_1, 2*pi);
        M_2 = mod(M_2, 2*pi);
    
        % Calculate the orbital period T
        T = 2 * pi * sqrt(a^3 / mu);
    
        % Calculate the mean anomaly difference
        delta_M = M_2 - M_1;
    
        % Calculate the time to travel between the true anomalies
        tof = (delta_M / (2 * pi)) * T;

    elseif e > 1
        sinhF_1 = (sqrt(e^2 - 1)*sin(nu_1)) / (1 + e*cos(nu_1));
        sinhF_2 = (sqrt(e^2 - 1)*sin(nu_2)) / (1 + e*cos(nu_2));

        F_1 = asinh(sinhF_1);
        F_2 = asinh(sinhF_2);

        tof = sqrt(-a^3 / mu)*((e*sinhF_2 - F_2) - (e*sinhF_1 - F_1));
    end
        
end