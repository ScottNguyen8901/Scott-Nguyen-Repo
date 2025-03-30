function gmst = GMST(date)
    %
    % DESCRIPTION
    %   Calculate the Greenwich Mean Sidereal Time (GMST) for a given UTC datetime.
    %
    % INPUTS    
    %   date_meas   (1,1)   datetime object []
    %
    % OUTPUTS
    %   gmst (1,1)   Greenwich Mean Sidereal Time (GMST) in radians []
    %
    % NOTES
    %   - The function uses a series approximation to compute GMST.
    %   - All angles, including GMST, are in radians.
    %
    % FUNCTION
    %
    % Define constants
    constants;

    % Extract components from datetime
    yy = year(date);
    mm = month(date);
    dd = day(date);
    hh = hour(date);
    min = minute(date);
    ss = second(date);

    % Calculate Julian Century from 2000
    century = (julian_date(yy, mm, dd, 0, 0, 0) - JD_2000) / CENTURY;

    % Convert time to total seconds
    total_seconds = hh * 3600 + min * 60 + ss;

    % Compute GMST using the series approximation
    gmst = ((-6.2E-6 * century + 0.093104) * century + 8640184.812866) * century + 24110.54841;
    
    % Adjust GMST for the time of day and mod 2*pi for range [0, 2*pi]
    gmst = mod(gmst * (pi / 43200) + w_E * total_seconds, 2 * pi);
end