function LMST_rad = calculate_LMST(year, month, day, hour, minute, longitude)
    % Function to calculate Local Mean Sidereal Time (LMST)
    % Input: year, month, day, hour (in MST), minute, longitude (in degrees)
    % Output: LMST in radians

    % Convert time to UT (MST is UTC-7)
    UT = hour + minute / 60 + 7;  % convert to UTC

    % Calculate Julian Date (JD)
    JD = 367 * year - floor(7 * (year + floor((month + 9) / 12)) / 4) ...
        + floor(275 * month / 9) + day + 1721013.5 + UT / 24;

    % Calculate T
    T = (JD - 2451545) / 36525;

    % Calculate Greenwich Mean Sidereal Time (GMST) in degrees
    GMST_deg = 280.46061837 + 360.98564736629 * (JD - 2451545) ...
        + T^2 * 0.000387933 - T^3 / 38710000;

    % Normalize GMST to [0, 360]
    GMST_deg = mod(GMST_deg, 360);

    % Convert GMST to radians
    GMST_rad = deg2rad(GMST_deg);

    % Convert longitude to radians
    longitude_rad = deg2rad(longitude);

    % Calculate Local Mean Sidereal Time (LMST) in radians
    LMST_rad = GMST_rad + longitude_rad;

    % Normalize LMST to [0, 2*pi]
    LMST_rad = mod(LMST_rad, 2 * pi);
    
    % Check if LMST is negative, convert to positive equivalent if necessary
    if LMST_rad < 0
        LMST_rad = LMST_rad + 2 * pi;
    end
end