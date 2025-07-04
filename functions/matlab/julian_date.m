function jd = julian_date(yy, mm, dd, hh, min, ss)
%
% DESCRIPTION
%   Calculate the Julian Date (JD) corresponding to the Gregorian date/time.
%
% INPUTS       size     Type       Description         Units
%   yy         (1,1)    Integer    Year (e.g., 2023)   []
%   mm         (1,1)    Integer    Month (1-12)        []
%   dd         (1,1)    Integer    Day of the month    []
%   hh         (1,1)    Integer    Hour (0-23)         []
%   min        (1,1)    Integer    Minute (0-59)       []
%   ss         (1,1)    Integer    Second (0-59)       []
%
% OUTPUTS      size     Type       Description         Units
%   jd         (1,1)    Double     Julian Date         [days]
%
% FUNCTION

    % Adjust month and year for calculation
    m = fix((mm - 14) / 12);
    
    % Calculate Julian Date
    jd = fix(dd - 32075 + (1461 * (yy + 4800 + m)) / 4) + ...
         fix(367 * (mm - 2 - 12 * m) / 12) - ...
         fix(3 * fix((yy + 4900 + m) / 100) / 4);
    
    % Adjust for time of day
    jd = jd + (hh - 12) / 24 + min / 1440 + ss / 86400;
end
