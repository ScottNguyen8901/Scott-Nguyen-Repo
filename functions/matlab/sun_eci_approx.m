function r_Sun_ECI_km = sun_eci_approx(date)
% Approximate Sun position in ECI (Earth-centered) in km.
% Good enough for day/night + eclipse gating; replace with higher-fidelity ephemeris if desired.

    % Julian Date
    jd = juliandate(date);

    % Days since J2000.0
    n = jd - 2451545.0;

    % Mean longitude (deg), mean anomaly (deg)
    L = mod(280.460 + 0.9856474*n, 360.0);
    g = mod(357.528 + 0.9856003*n, 360.0);

    % Ecliptic longitude (deg)
    lambda = L + 1.915*sind(g) + 0.020*sind(2*g);

    % Obliquity of ecliptic (deg)
    eps = 23.439 - 0.0000004*n;

    % Distance in AU (approx)
    R = 1.00014 - 0.01671*cosd(g) - 0.00014*cosd(2*g);

    % Convert to ECI (assume mean equator of date approx)
    x = R*cosd(lambda);
    y = R*cosd(eps)*sind(lambda);
    z = R*sind(eps)*sind(lambda);

    AU_km = 149597870.7; % [km]
    r_Sun_ECI_km = AU_km * [x; y; z];
end