function theta_g = calc_GST(date)
%
% DESCRIPTION
%   Calculate Greenwich Sidereal Time (GST) for a given date (in datetime format)
%   and return it in radians.
%
% INPUTS     size     Type       Description                         Units
%   date     (1,1)    datetime   Date and time for the calculation   []
%
% OUTPUTS    size     Type       Description                     Units
%   theta_g  (1,1)    Double     Greenwich Sidereal Time         [rad]
%
% NOTES
%
% FUNCTION
    
    % Convert the datetime to Julian Date (JD)
    Y = year(date);
    M = month(date);
    D = day(date);
    UT = hour(date) + minute(date)/60 + second(date)/3600;
    
    % Julian Date formula for the given UTC time
    JD = 367 * Y - floor(7 * (Y + floor((M + 9) / 12)) / 4) + ...
         floor(275 * M / 9) + D + 1721013.5 + UT / 24;
    
    % Calculate Julian centuries from J2000.0
    T = (JD - 2451545.0) / 36525;
    
    % Calculate the Greenwich Mean Sidereal Time (GMST) in degrees
    GMST_deg = 280.46061837 + 360.98564736629 * (JD - 2451545) + T^2 * 0.000387933 - T^3 / 38710000;
    theta_g = mod(GMST_deg, 360); % Ensure GMST is within 0 to 360 degrees
    
end