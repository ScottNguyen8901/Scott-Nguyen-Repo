function f = anom_to_true_anom(anomaly, e)
    %
    % DESCRIPTION
    %   Converts the given orbital anomaly (eccentric, hyperbolic, or parabolic) 
    %   to the true anomaly for a given eccentricity.
    %
    % INPUTS         size         Type     Description                               Units
    %   anomaly      (1,1)        Double   Orbital anomaly (eccentric, hyperbolic, or parabolic) [rad]
    %   e             (1,1)        Double   Eccentricity of the orbit                [unitless]
    %
    % OUTPUTS        size         Type     Description                               Units
    %   theta        (1,1)        Double   True anomaly of the orbit                [rad]
    %
    % NOTES
    %   This function calculates the true anomaly based on the provided orbital 
    %   anomaly and eccentricity. It applies different formulas depending on the type 
    %   of anomaly (eccentric, hyperbolic, or parabolic).
    %

    if e < 1  % Elliptical Orbit
        % Convert eccentric anomaly to true anomaly for elliptical orbit
        f = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(anomaly / 2));
        
    elseif e > 1  % Hyperbolic Orbit
        % Convert hyperbolic anomaly to true anomaly for hyperbolic orbit
        f = 2 * atan(sqrt((e + 1) / (e - 1)) * tanh(anomaly / 2));
        
    else  % Parabolic Orbit (e = 1)
        % Convert parabolic anomaly to true anomaly for parabolic orbit
        f = 2 * atan(sqrt(2) * anomaly);
    end
end