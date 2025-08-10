function a = ground_track_to_sma(N1, N2, mu)
    %
    % DESCRIPTION
    %   Calculate the semi-major axis of an orbit based on the difference in the 
    %   longitudes of ascending nodes and the gravitational parameter.
    %
    % INPUTS    size    Type
    %   N1      (1,1)   Double  Longitude of Ascending Node 1 [deg]
    %   N2      (1,1)   Double  Longitude of Ascending Node 2 [deg]
    %   mu      (1,1)   Double  Gravitational parameter       [km^3/s^2]
    %
    % OUTPUTS   size    Type
    %   a       (1,1)   Double  Semi-major axis of orbit [km]
    %
    % FUNCTION 

    % Calculate the difference in longitudes
    delta_N = mod(N2 - N1 + 360, 360);
    
    % Calculate the time period in hours
    TP_hours = 24 * (1 - delta_N / 360);
    
    % Convert time period to seconds
    TP_seconds = 3600 * TP_hours;

    % Calculate the semi-major axis
    a = (mu * (TP_seconds / (2 * pi))^2)^(1/3);
    
end
