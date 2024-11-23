clear
clc

% Dimensions and masses
a = 1; % Prism base width
b = 1; % Prism base height
h = 2; % Prism height
M = 100; % Prism mass

L = 2; % Panel length
w = 0.6; % Panel width
m = 10; % Panel mass

% Offset for panels
r = (a + L) / 2;

% Rectangular prism moment of inertia
I_rp = (M / 12) * diag([b^2 + h^2, a^2 + h^2, a^2 + b^2]);

% Panel moment of inertia
I_p = (m / 12) * diag([w^2 + L^2, w^2, w^2]);

% Offset vector for first panel
r_vec1 = [0; -r; 0];
I_panel1 = I_p + m * (r^2 * eye(3) - r_vec1 * r_vec1');

% Offset vector for second panel
r_vec2 = [0; r; 0];
I_panel2 = I_p + m * (r^2 * eye(3) - r_vec2 * r_vec2');

% Total moment of inertia tensor
I_o = I_rp + I_panel1 + I_panel2;

%%
clear
clc

I = [51.6294, 1.8192, 0.4055;
     1.8192, 41.2647, 0.6444;
     0.4055, 0.6444, 31.9150];

[V,D]=eig(I)

%%
% Given moment of inertia matrix
I = [181.2123, 9.1086, 10.6243;
     9.1086, 192.0960, 11.1512;
     10.6243, 11.1512, 44.8092];

% Calculate the eigenvalues of the matrix
eigenvalues = eig(I);

% The minor axis moment of inertia is the smallest eigenvalue
minor_axis_moi = min(eigenvalues);

% Display the result
disp(['Minor axis moment of inertia: ', num2str(minor_axis_moi)]);

%%
% Given moment of inertia matrix
I = [181.2123, 9.1086, 10.6243;
     9.1086, 192.0960, 11.1512;
     10.6243, 11.1512, 44.8092];

% Calculate the eigenvalues of the matrix
eigenvalues = eig(I);

% The major axis moment of inertia is the largest eigenvalue
major_axis_moi = max(eigenvalues);

% Display the result
disp(['Major axis moment of inertia: ', num2str(major_axis_moi)]);

%%
% Given moment of inertia matrix
I = [181.2123, 9.1086, 10.6243;
     9.1086, 192.0960, 11.1512;
     10.6243, 11.1512, 44.8092];

% Calculate the eigenvalues and eigenvectors of the matrix
[eigenvectors, eigenvalues] = eig(I);

% Find the index of the smallest eigenvalue (minor axis)
[~, minor_index] = min(diag(eigenvalues));

% The corresponding eigenvector is the unstable spin axis
unstable_spin_axis = eigenvectors(:, minor_index);

% Display the result
disp('Unstable spin axis:');
disp(unstable_spin_axis);

%%

clear
clc

% Dimensions and masses
a = 1; % Prism base width
b = 1; % Prism base height
h = 2; % Prism height
M = 100; % Prism mass

L = 2; % Panel length
w = 0.6; % Panel width
m = 10; % Panel mass

% Offset for panels
r = (a + L) / 2;

% Rectangular prism moment of inertia
I_rp = (M / 12) * diag([b^2 + h^2, a^2 + h^2, a^2 + b^2]);

% Panel moment of inertia
I_p = (m / 12) * diag([w^2 + L^2, w^2, w^2]);

% Offset vector for first panel
r_vec1 = [0; -r; 0];
I_panel1 = I_p + m * (r^2 * eye(3) - r_vec1 * r_vec1');

% Offset vector for second panel
r_vec2 = [0; r; 0];
I_panel2 = I_p + m * (r^2 * eye(3) - r_vec2 * r_vec2');

% Total moment of inertia tensor
I_o = I_rp + I_panel1 + I_panel2;