close
clear
clc

% Given quantities
% Cylinder
d = 1;            % [m]
r = d / 2;        % [m]
h = 2;            % [m]
mu = 100;         % [kg / m^3]
V = pi * r^2 * h; % [m^3]
m_cyl = mu * V;   % [kg]

% Rod
l = 2;     % [m]
m_rod = 1; % [kg]

% Panel
L = 2;               % [m]
w = 0.5;             % [m]
m_pan = 4;           % [kg]
theta = deg2rad(87); % [rad]

% Mass centers of components about center of cylinder
r_cyl_c_b = [0; 0; 0];
r_rod_c_b = [0; r + l / 2; h];
r_pan_c_b = [0; r + l / 2; h];

% Moment of Inertia Calculations about mass centers
I_cyl_c_b = m_cyl * diag([(3 * r^2 + h^2) / 12; ...
                      (3 * r^2 + h^2) / 12; ...
                      r^2 / 2]);

I_rod_c_b = m_rod * diag([l^2 / 12, 0, l^2 / 12]);

I_pan_c_a = m_pan * diag([(L + w)^2 / 12; ...
                          w^2 / 12; ...
                          L^2 / 12]);

% Change of Base Vector
R_ba = rot_2(theta);
I_pan_c_b = R_ba * I_pan_c_a * R_ba';

% Parallel Axis Theorem
I_cyl_o_b = I_cyl_c_b - m_cyl * skew(r_cyl_c_b) * skew(r_cyl_c_b);
I_rod_o_b = I_rod_c_b - m_rod * skew(r_rod_c_b) * skew(r_rod_c_b);
I_pan_o_b = I_pan_c_b - m_pan * skew(r_pan_c_b) * skew(r_pan_c_b);

% Total Moment of Inertia about Center of Cylinder
I_comp_o_b = I_cyl_o_b + I_rod_o_b + I_pan_o_b;

% Location of center of mass about point o 
M = m_cyl + m_rod + m_pan;
r_c_o = (m_cyl * r_cyl_c_b + m_rod * r_rod_c_b + m_pan * r_pan_c_b) / M;

% Total Moment of Inertia about Center of Mass
I_comp_c_b = I_comp_o_b - M * skew(r_c_o) * skew(r_c_o);

[axes_new, mois_new] = reorder_axes_and_mois(I_comp_c_b);

% Converting rotation matrix to euler angle and axis
[phi, a] = rot_to_eul(axes_new);

phi = rad2deg(phi);