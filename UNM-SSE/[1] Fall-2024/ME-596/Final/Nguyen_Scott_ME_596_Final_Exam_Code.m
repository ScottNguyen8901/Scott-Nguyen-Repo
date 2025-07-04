clear
close
clc

R = 6378;
mu = 3.986e5; 

%% Problem 1

% A gravity-gradient stabilized satellite uses its attitude determination 
% sensors to determine that the nadir vector, o3, has components, 
% before normalization,  (0   0.207   0.802) in Fb.  What is the angle 
% between the b3 direction and the direction to the Earth?  
% Give your answer in radians with four digits to the right of the decimal place.

clear

o_3 = [0; 0.207; 0.802];
o_3 = o_3 / norm(o_3);
b_3 = [0; 0; 1];

theta = acos(dot(o_3, b_3) / (norm(o_3) * norm(b_3)));

fprintf('Problem 1: %0.4f rad\n', theta)


%% Problem 2

% A torque-free axisymmetric body has principal moments of inertia satisfying 
% C=3A/4.  At t=0, the precession and spin angles are both zero, and the 
% angular velocity vector expressed in the body frame has components (1,2,19.5) rad/s.
% What is the nutation angle?  Give your answer in radians with four digits
%  to the right of the decimal place.

clear

clear

% Define symbolic variable for A
syms A;

% Given angular velocity vector (rad/s)
omega = [1; 2; 19.5];

% Inertia matrix for axisymmetric body
C = (3/4) * A;  % C = 3A/4
I = [A, 0, 0; 
     0, A, 0; 
     0, 0, C];  % Inertia matrix

% Angular momentum vector (L = I * omega)
L = I * omega;

% Compute the magnitude of the angular momentum vector
L_magnitude = sqrt(L(1)^2 + L(2)^2 + L(3)^2);

% Compute the z-component of the angular momentum vector (Lz)
Lz = L(3);

% Calculate the nutation angle (theta)
cos_theta = Lz / L_magnitude;
theta = acos(cos_theta);  % Nutation angle in radians

% Simplify the expression for theta
theta_simplified = simplify(theta, 'IgnoreAnalyticConstraints', true);

% Convert to numerical value (substitute A = 1 or any specific value for A)
theta_numeric = double(subs(theta_simplified, A, 1));  % Example for A=1

% Display the simplified result using fprintf with 4 decimal places
fprintf('\nProblem 2: Nutation Angle %.4f radians\n', theta_numeric);

%% Problem 3

% A torque-free axisymmetric body has principal moments of inertia 
% satisfying C=3A/4.  At t=0, the angular velocity vector expressed in the 
% body frame has components (-1,-4,20) rad/s. What is the phase angle in 
% the solution for ω1 and ω2?  Give your answer in radians, 
% between 0 and 2π, with four digits to the right of the decimal place.

clear

w = [-1; -4; 20];
phi = mod(atan2(w(1), w(2)), 2 * pi);
fprintf('\nProblem 3: Phase angle %.4f rad.\n', phi);


%% Problem 4

% A torque-free rigid body is initially rotating with angular velocity 
% (0.9572   0.4854   0.974) rad/s, expressed in the principal Fb.  
% The principal moments of inertia in this frame are 200, 150, and 90 kg m2.
% The rigid body is subject to internal energy dissipation that does 
% not discernibly affect the moments of inertia.  What angular velocity 
% will the body eventually reach as t goes to infinity?  
% Give your answer as a scalar, in rad/s, with four digits to the 
% right of the decimal place.

clear

w0 = [0.9572; 0.4854; 0.974];  % Initial angular velocity (rad/s)
I = diag([200; 150; 90]);     % Moments of inertia (kg*m^2)


% Define the ODE system for angular velocity
dw_dt = @(t, w) -inv(I) * skew(w) * I * w;

% Set the time span for the integration (large t to approach infinity)
t_span = [0, 1e5];  % A large time span to simulate long-term behavior

% Integrate the system using ode45
[t, w_t] = ode45(@(t, w) dw_dt(t, w), t_span, w0);

% The final angular velocity at large t will be the last value of w_t
final_w = norm(w_t(end, :));  % Compute the scalar magnitude of the final angular velocity

% Display the result
fprintf('\nProblem 4:  %.4f rad/s\n', final_w);

w_norm = vecnorm(w_t, 2, 2);
figure;
plot(t, w_norm, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Norm of Angular velocity (rad/s)');
title('Norm of Angular velocity over time');
grid on;

%% Problem 5

% If you use the Triad algorithm to estimate the rotation matrix Rbi, 
% using the following data, and using v1 as the "exact" measurement, 
% then what is the value of Rbiv1i?


clear

v_1i = [0.8739; 0.1549; 0.4606];  
v_2i = [0.5927; 0.5127; 0.6210]; 

v_1b = [0.6763; 0.4915; 0.5487];
v_2b = [0.6973; 0.3680; 0.6151];  

v_1i = v_1i / norm(v_1i);
v_2i = v_2i / norm(v_2i);

v_1b = v_1b / norm(v_1b);
v_2b = v_2b / norm(v_2b);

v_b = [v_1b, v_2b];
v_i = [v_1i, v_2i];
choice = 1;

[R_bi_tri, J] = triad(v_b, v_i, choice);

ans = R_bi_tri * v_1i;

fprintf('\nProblem 5: Rbi * v1i = [%0.4f %0.4f %0.4f]^T\n', ans(1), ans(2), ans(3));

%% Problem 6

% The moment of inertia matrix for a spacecraft about its mass center 
% in body-fixed reference frame Fb is: Find the principal reference frame 
% such that the frame's unit vectors are as close as possible to 
% alignment with  Fb's unit vectors.  What is the angle between the b1 
% vector and the principal axis unit vector that is closest to the b1 
% vector?  Give your answer in radians with 4 digits to the right of 
% the decimal place

clear

I = [99, -3, -3;...
    -3, 118, 5;...
    -3, 5, 91];

[axes_new, mois_new] = reorder_axes_and_mois(I);

b_1 = [1; 0; 0];
b_1p = axes_new(:,1);

theta = acos(dot(b_1, b_1p) / (norm(b_1) * norm(b_1p)));

fprintf('\nProblem 6: %0.4f rad\n', theta)

%% Problem 7

% The moment of inertia matrix for a 20-kg spacecraft about a specific 
% point a in body-fixed reference frame Fb is (in kg m2) The position 
% vector from the mass center o to the point a has components in Fb of 
% [0.2  0  0.3], with units of meters. What is the largest principal moment 
% of inertia about the mass center?  Give your answer in kg m2, with four 
% digits to the right of the decimal place.

clear

% Given data
m = 20;
I_a = [138, -2, -7; 
       -4, 100, -8;
       -7, -7, 105];

r_oa = [0.2; 0; 0.3];

I_o = I_a  m * skew(r_oa) * skew(r_oa);
[axes, mois] = eig(I_o);
moi_max = max(diag(mois));

fprintf('\nProblem 7: %0.4f kg * m^2\n', moi_max)