function nu = ecc_to_true_anom(E, e)
    %
    % DESCRIPTION
    %   Convert Eccentric Anomaly to True Anomaly for an elliptical orbit.
    %   The function takes the Eccentric Anomaly (E) and orbit eccentricity (e)
    %   as inputs and computes the corresponding True Anomaly (nu).
    %
    % INPUTS         size         Type    Description                   Units
    %   E            (1,1)        Double  Eccentric Anomaly             [rad]
    %   e            (1,1)        Double  Eccentricity of the orbit     []
    %
    % OUTPUTS        size         Type    Description                   Units
    %   nu           (1,1)        Double  True Anomaly                  [rad]
    %
    % NOTES
    %   The True Anomaly (nu) is the angle between the direction of periapsis
    %   and the current position of the orbiting body, measured at the focus
    %   of the ellipse.
    %
    % FUNCTION

    % Calculate True Anomaly using the relationship between Eccentric Anomaly and True Anomaly
    nu = 2 * atan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2));

end
