clear
close
clc

%% Problem 1

% An axisymmetric spacecraft has principal moments of inertia (50, 50, 70) kg m2.  
% At t=0, the spacecraft is rotating with angular velocity components (1,1,10) rad/s.   
% What is the value of the constant c1 in the solution to the linear differential equation for w1?
% Give your answer in rad/s, with four digits to the right of the decimal place.

% Given parameters
I = [50; 50; 70]; % Moments of inertia
t = 0; % Initial time
w = [1, 1, 10]; % Angular velocity components

% Constants
A = I(1);
C = I(3);

W = w(3); % Angular velocity about the z-axis
W_hat = W * (A - C) / A; % Nutation frequency

% Define symbolic variables
syms c_1 c_2

% Define the equation
eqn = w(1) == c_1 * cos(W_hat * t) + c_2 * sin(W_hat * t);

% Solve for c_1 at t = 0
c_1_solution = solve(eqn, c_1);

% Display the result
fprintf('Problem 1: %.4f\n', double(c_1_solution));

%% Problem 2

% An asymmetric rigid body is rotating without any applied torque, and has 
% a slight wobble about its nominal spin axis.  All numbers are given in 
% metric units (kg, m, s). The angular momentum vector has the following 
% components in the body-fixed principal reference frame:  [50  0  700].
% The rotational kinetic energy is 3525, and it is known that I3 = 70.
% What is the "1" component of the angular velocity vector in the body frame?  
% Give your answer with 4 digits to the right of the decimal place

% Given data
H = [50, 0, 700];  % Angular momentum components [H1, H2, H3]
T = 3525;          % Rotational kinetic energy
I3 = 70;           % Moment of inertia about axis 3

% Extract angular momentum components
H1 = H(1);
H2 = H(2);
H3 = H(3);

% Calculate omega_3
omega3 = H3 / I3;

% Symbolically solve for I1, I2, omega1, and omega2
syms I1 I2 omega1 omega2
eq1 = 0.5 * (I1 * omega1^2 + I2 * omega2^2 + I3 * omega3^2) == T; % Energy equation
eq2 = H1 == I1 * omega1;                                         % Angular momentum for axis 1
eq3 = H2 == I2 * omega2;                                         % Angular momentum for axis 2 (H2 = 0)

% Solve equations
solution = solve([eq1, eq2, eq3], [I1, I2, omega1, omega2], 'Real', true);

% Extract solutions
I1_val = double(solution.I1);
I2_val = double(solution.I2);
omega1_val = double(solution.omega1);
omega2_val = double(solution.omega2);

% Display results
fprintf('\nProblem 2:\n');
fprintf('I1 = %.4f\nI2 = %.4f\nomega_1 = %.4f rad/s\nomega_2 = %.4f rad/s\nomega_3 = %.4f rad/s\n', ...
        I1_val, I2_val, omega1_val, omega2_val, omega3);

%% Problem 3

% A vector has components in Fb of (5 8 2).  What is the angle between 
% the vector and the b2 axis? Give your answer in radians with four digits
% to the right of the decimal place.

v = [5; 8; 2];
b_2 = [0; 1; 0];

theta = acos(dot(v, b_2) / (norm(v) * norm(b_2)));

fprintf('\nProblem 3: %.6f rad\n', theta)

%% Problem 4

% A spacecraft has moment of inertia matrix in Fa with elements as follows: 
% What is the principal moment of inertia about the "2" axis 
% (where the principal reference frame is chosen so that the axes are 
% "close" to the corresponding axes of the Fa)?  Give your answer with four 
% digits to the right of the decimal place.

I = [50.0124, -0.0471, 0.4966;
    -0.0471, 50.1781, -1.8783;
     0.4966, -1.8783, 69.8095];

[axes_new, mois_new] = reorder_axes_and_mois(I);

fprintf('\nProblem 4: I_2 %.4f\n', mois_new(2,2))

%% Problem 5

% A low-accuracy attitude determination system has determined that
% mi = [ 0.2 0.8 0.4 ]T, mb = [ 0.6 0.6 0.3 ]T, si = [0.5 0.7 0.9 ]T, and sb = [0.4 0.7 0.5 ]T
% Use the q-method to determine the maximum eigenvalue of the K matrix, using w1=1/2 and w2=1/2.   
% Do not forget to normalize the unit vectors!
% Enter your value of λmax with four digits to the right of the decimal place.

m_i = [0.2; 0.8; 0.4];
m_b = [0.6; 0.6; 0.3];

s_i = [0.5; 0.7; 0.9];
s_b = [0.4; 0.7; 0.5];

m_i = m_i / norm(m_i);
m_b = m_b / norm(m_b);

s_i = s_i / norm(s_i);
s_b = s_b / norm(s_b);

v_b = [m_b, s_b];
v_i = [m_i, s_i];
w_k = [1/2, 1/2];

[q, R_bi_q, J, eig_max] = q_method(v_b, v_i, w_k);

fprintf('\nProblem 5: Max eigen-value %.12f\n', eig_max)

%% Problem 6

% A given rigid body has the following moment of inertia matrix expressed in Fa:
% What is the largest of the principal moments of inertia?  Give your 
% answer with four digits to the right of the decimal place.

I = [51.6294, 1.8192, 0.4055;
     1.8192, 41.2647, 0.6444; 
     0.4055, 0.6444, 31.9150];

[axes_new, mois_new] = reorder_axes_and_mois(I);

fprintf('\nProblem 6: Max MOI %.4f\n', max(mois_new(:)))