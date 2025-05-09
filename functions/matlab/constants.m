% DESCRIPTION
%   Define physical and astronomical constants used in orbital mechanics
%   and Earth rotation calculations.
%
% CONSTANTS      Value                  Description                            Units
%   R_E          6378.1370              Earth's equatorial radius              [km]
%   mu_E         3.986004418E5          Earth's gravitational parameter        [km^3/s^2]
%   w_E          7.2921150E-5           Earth's rotational rate (scalar)       [rad/s]
%   omega_E_vec  [0; 0; 7.2921150e-5]   Earth's rotational vector              [rad/s]
%   CENTURY      36525                  Days in a Julian century               [days]
%   JD_2000      2451545.0              Julian Date of J2000 epoch             []
%   f_E          1/298.257223563        Earth's flattening factor              []
%   sec_day      86400                  Number of seconds in a day             [s]
%   J_2          0.00108263             Earth's J2 coefficient                 []
%   Cd           2.2                    Drag coefficient                       []
%   A            20                     Area normal to velocity vector         [m^2]
%   m            400                    Spacecraft mass                        [kg]
%   r_s_TH       [0;0;R_E]              Position of sensor in TH frame         [km]

R_E         = 6378.1370;           % km
mu_E        = 3.986004418E5;       % km^3/s^2
w_E         = 7.2921150E-5;        % rad/s
omega_E_vec = [0; 0; w_E];         % rad/s
CENTURY     = 36525;               % days
JD_2000     = 2451545.0;           % JD of 01 January, 2000 12h UTC
f_E         = 1/298.257223563;     % []
sec_day     = 86400;               % s
J_2         = 0.00108263;          % []

Cd          = 2.2;                 % []
A           = 0.03;                % m^2
m           = 5;                   % kg

r_s_TH      = [0; 0; R_E];         % [km]
