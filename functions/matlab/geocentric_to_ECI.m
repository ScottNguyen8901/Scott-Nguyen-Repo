function [x_ECI, theta] = geocentric_to_ECI(r_TH, lat, lon, date)
%
% DESCRIPTION
%   Converts the sensor's geocentric position and velocity to Earth-Centered Inertial (ECI) 
%   coordinates and calculates the Greenwich Mean Sidereal Time (GMST) for a given location 
%   and time.
%
% INPUTS         size     Type       Description                      Units
%   lat          (1,1)    Double     Sensor latitude (North positive) [rad]
%   lon          (1,1)    Double     Sensor longitude (East positive) [rad]
%   date_0       (1,1)    datetime   Date and time                    []
%
% OUTPUTS        size     Type       Description                      Units
%   r_sensor_ECI (9,1)    Double     Sensor ECI state [pos; vel; acc] [DU, DU/TU, DU/TU^2]
%   theta        (1,1)    Double     Local mean sidereal time         [rad]
%
% NOTES
%
% FUNCTION

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
    p = R_E * cos(lat);
    ct = cos(theta);
    st = sin(theta);
    
    % Velocity of the sensor in ECI coordinates
    x_ECI(4:6,1) = [-st; ct; 0] * p * w_E;
    
    % Acceleration of the sensor in ECI coordinates
    x_ECI(7:9,1) = [-ct; -st; 0] * p * w_E^2;
    end