function [x_ECI, theta] = geocentric_to_ECI(r_TH, lat, lon, date)
    % GEOCENTRIC_TO_ECI Converts the sensor's geocentric position and velocity
    % to Earth-centered inertial (ECI) coordinates and calculates the Greenwich 
    % Mean Sidereal Time (GMST) for a given location at a specified time.
    %
    % Inputs:
    %   lat    - Latitude of the sensor in radians (North positive)
    %   lon    - Longitude of the sensor in radians (East positive)
    %   date_0 - datetime object specifying the date and time
    %
    % Outputs:
    %   r_sensor_ECI - 9x1 state vector containing the sensor's position, velocity, 
    %                  and acceleration in ECI coordinates [position; velocity; acceleration]
    %   theta        - Local mean sidereal time (LST) in radians
    %
    % Notes:

    % Constants
    constants;
    
    % Calculate Greenwich Mean Sidereal Time (GMST)
    theta_g = mod(GMST(date), 2*pi);
    
    % Local Sidereal Time (LST) is the longitude plus GMST
    theta = lon + theta_g;
    
    % Rotation matrix to convert from TH coordinates to ECI coordinates
    Q_TH_ECI = rot_ijk_sez(lat, theta);
    
    % Convert sensor's position from TH to ECI coordinates
    x_ECI = Q_TH_ECI' * r_TH;
    
    % Calculate the sensor's velocity and acceleration in ECI coordinates
    p = R_E * cos(lat);  % Projection of the sensor's distance in the ECI plane
    ct = cos(theta);   % Cosine of the local sidereal time
    st = sin(theta);   % Sine of the local sidereal time
    
    % Velocity of the sensor in ECI coordinates
    x_ECI(4:6,1) = [-st; ct; 0] * p * w_E;
    
    % Acceleration of the sensor in ECI coordinates
    x_ECI(7:9,1) = [-ct; -st; 0] * p * w_E^2;
    end