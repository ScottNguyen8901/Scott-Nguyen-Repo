function tof = tof_ta(a, e, nu_1, nu_2, mu)
    %
    % DESCRIPTION
    %   Compute the time of flight between two true anomalies for an elliptical orbit.
    %   The function uses the eccentricity, semi-major axis, gravitational parameter, and
    %   the true anomalies (in degrees) to compute the time taken to travel between them.
    %
    % INPUTS    size     Type    Description                   Units   
    %   e       (N,1)    Double  Eccentricity                  []
    %   a       (N,1)    Double  Semi-major axis               [DU]
    %   mu      (1,1)    Double  Gravitational parameter       [DU^3/TU^2]
    %   nu_1    (N,1)    Double  True anomaly 1                [rad]
    %   nu_2    (N,1)    Double  True anomaly 2                [rad]
    %
    % OUTPUTS   size     Type    Description                   Units
    %   tof     (N,1)    Double  TOF between the two TA        [TU]
    %
    % NOTES
    %

    % Initialize TOF array
    tof = zeros(size(a));

    % Loop through each element to calculate TOF
    for i = 1:length(a)
        if e(i) > 0 && e(i) < 1
            % Calculate the eccentric anomaly E for both true anomalies
            tan_nu_1_over_2 = tan(nu_1(i) / 2);
            sqrt_term = sqrt((1 + e(i)) / (1 - e(i)));
            tan_E1_over_2 = tan_nu_1_over_2 / sqrt_term;
            E_1 = 2 * atan(tan_E1_over_2);

            tan_nu_2_over_2 = tan(nu_2(i) / 2);
            tan_E_2_over_2 = tan_nu_2_over_2 / sqrt_term;
            E_2 = 2 * atan(tan_E_2_over_2);

            % Calculate the mean anomalies M from E
            M_1 = E_1 - e(i) * sin(E_1);
            M_2 = E_2 - e(i) * sin(E_2);

            % Normalize the mean anomalies M within the range [0, 2*pi]
            M_1 = mod(M_1, 2*pi);
            M_2 = mod(M_2, 2*pi);

            % Calculate the orbital period T
            T = 2 * pi * sqrt(a(i)^3 / mu);

            % Calculate the mean anomaly difference
            delta_M = M_2 - M_1;

            % Calculate the time to travel between the true anomalies
            tof(i) = (delta_M / (2 * pi)) * T;

        elseif e(i) > 1
            sinhF_1 = (sqrt(e(i)^2 - 1) * sin(nu_1(i))) / (1 + e(i) * cos(nu_1(i)));
            sinhF_2 = (sqrt(e(i)^2 - 1) * sin(nu_2(i))) / (1 + e(i) * cos(nu_2(i)));

            F_1 = asinh(sinhF_1);
            F_2 = asinh(sinhF_2);

            tof(i) = sqrt(-a(i)^3 / mu) * ((e(i) * sinhF_2 - F_2) - (e(i) * sinhF_1 - F_1));
        end
    end
end
