function nu = ecc_to_true_anom(E, e)
    % DESCRIPTION
    %   Convert Eccentric Anomaly (E) or Hyperbolic Anomaly (H) to True Anomaly (nu).
    %   The function takes the Anomaly (E or H) and orbit eccentricity (e)
    %   as inputs and computes the corresponding True Anomaly (nu).
    %
    % INPUTS         size         Type    Description                   Units
    %   E            (1,1)        Double  Eccentric or Hyperbolic Anomaly [rad]
    %   e            (1,1)        Double  Eccentricity of the orbit       []
    %
    % OUTPUTS        size         Type    Description                   Units
    %   nu           (1,1)        Double  True Anomaly                  [rad]
    %
    % NOTES
    %   For elliptical orbits (0 <= e < 1), E represents the Eccentric Anomaly.
    %   For hyperbolic orbits (e > 1), E represents the Hyperbolic Anomaly (H).
    %   The True Anomaly (nu) is the angle between the direction of periapsis
    %   and the current position of the orbiting body, measured at the focus
    %   of the conic section.
    %
    % FUNCTION

    if e < 1
        % Elliptical orbit: Calculate True Anomaly using Eccentric Anomaly
        nu = 2 * atan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2));
    elseif e > 1
        % Hyperbolic orbit: Calculate True Anomaly using Hyperbolic Anomaly
        nu = 2 * atan2(sqrt(e + 1) * sinh(E / 2), sqrt(e - 1) * cosh(E / 2));
    else
        error('Parabolic orbits (e = 1) are not supported by this function.');
    end
end
