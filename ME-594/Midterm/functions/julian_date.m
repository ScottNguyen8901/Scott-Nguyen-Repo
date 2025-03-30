function jd = julian_date(yy, mm, dd, hh, min, ss)
    %
    % DESCRIPTION
    %   Calculate the Julian Date (JD) corresponding to the Gregorian date/time.
    %
    % INPUTS    
    %   yy    (1,1)   Year (e.g., 2023) [integer]    
    %   mm    (1,1)   Month (1-12) [integer]
    %   dd    (1,1)   Day of the month [integer]
    %   hh    (1,1)   Hour (0-23) [integer]
    %   min   (1,1)   Minute (0-59) [integer]
    %   ss    (1,1)   Second (0-59) [integer]
    %
    % OUTPUTS
    %   jd    (1,1)   Julian Date corresponding to the input date/time [days]
    %
    % NOTES
    %   The calculation assumes the Gregorian calendar is in use.
    %   The Julian Date is returned as a decimal value representing the number of
    %   days since the Julian epoch (January 1, 4713 BCE, 12:00 UT).
    %
    % FUNCTION
    %

    % Adjust month and year for calculation
    m = fix((mm - 14) / 12);
    
    % Calculate Julian Date
    jd = fix(dd - 32075 + (1461 * (yy + 4800 + m)) / 4) + ...
         fix(367 * (mm - 2 - 12 * m) / 12) - ...
         fix(3 * fix((yy + 4900 + m) / 100) / 4);
    
    % Adjust for time of day
    jd = jd + (hh - 12) / 24 + min / 1440 + ss / 86400;
end
