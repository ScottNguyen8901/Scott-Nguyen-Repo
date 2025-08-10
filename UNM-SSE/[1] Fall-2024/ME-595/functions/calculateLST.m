function LST = calculateLST(dateVec, utcOffset, longitude)
    % Calculate Local Sidereal Time (LST)
    % Input:
    %   dateVec   - [year, month, day, hour, minute, second]
    %   utcOffset - Time difference from UTC in hours
    %   longitude - Longitude of the location in degrees (+ for East, - for West)
    % Output:
    %   LST       - Local Sidereal Time in degrees
    
    % Convert date vector to Julian Date
    Y = dateVec(1);
    M = dateVec(2);
    D = dateVec(3);
    UT = dateVec(4) + dateVec(5) / 60 + dateVec(6) / 3600; % Convert to decimal hours
    
    if M <= 2
        Y = Y - 1;
        M = M + 12;
    end
    
    % Julian Date formula for the given UTC time
    JD = 367 * Y - floor(7 * (Y + floor((M + 9) / 12)) / 4) + ...
         floor(275 * M / 9) + D + 1721013.5 + UT / 24;
    
    % Calculate Julian centuries from J2000.0
    T = (JD - 2451545.0) / 36525;
    
    % Updated GMST formula for high accuracy (IAU 2000 model)
    GMST = 280.46061837 + 360.98564736629 * (JD - 2451545.0) + ...
           T^2 * 0.000387933 - T^3 / 38710000;
    
    % Normalize GMST to the range [0, 360) degrees
    GMST = mod(GMST, 360);
    
    % Calculate LST by adjusting for the observer's longitude and UTC offset
    LST = GMST + longitude + (utcOffset * 15);
    
    % Normalize LST to the range [0, 360) degrees
    LST = mod(LST, 360);
    
    % Ensure LST is positive (normalize the angle)
    if LST < 0
        LST = LST + 360;
    end
end
