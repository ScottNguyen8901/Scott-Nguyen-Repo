clear
close
clc

%% Generate SS 

% Define the system matrices for a highly coupled system
A = [-2500, 2;   % System matrix with stronger coupling, adjusting to match output ranges
     50, -4];   % The system is designed to produce outputs in the desired range

B = [5, 2;   % Adjusted Input matrix for scaling
     5, 500];     % The input matrix has been adjusted for better scaling

C = [1, 0;     % Output matrix to observe both states
     0, 1];

D = [0, 0;     % No direct feedthrough
     0, 0];

% Define the system
ss_mat = ss(A, B, C, D);

%% HESO
% Define the parameters
m = 3;
B0 = [5, 0; 0, 100];
U = [50; 60];
w_o = 10;

% Call the gen_HESO function
[A_bar, B_bar, Beta] = gen_HESO(m, B0, U, w_o);

model = 'test';
simout = sim(model);

%% Plotting

y = simout.Y.signals.values;
t = simout.Y.time;

MR = y(:, 1);
thrust = y(:, 2);

figure;

subplot(1, 2, 1);
plot(t, MR, 'r', 'LineWidth', 2);
title('Mixture Ratio');
xlabel('Time [s]');
ylabel('Mixture Ratio');
axis square;
grid on;

subplot(1, 2, 2);
plot(t, thrust, 'b', 'LineWidth', 2);
title('Thrust');
xlabel('Time [s]');
ylabel('Thrust');
axis square;
grid on;