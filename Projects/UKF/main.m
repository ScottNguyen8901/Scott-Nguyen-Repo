% demo_UKF_DMC.m
% UKF demos: (1) noise only, (2) noise + sinusoidal perturbation,
% and (3) DMC-augmented UKF tracking the perturbation as a Singer process.

%% Setup and initial configuration
close all; clear; clc;
rng(1);   % Reproducibility

% --- Plotting properties ---
fs = 12; fn = 'Times New Roman'; lw = 2;
set(groot, 'DefaultAxesFontSize', fs, 'DefaultAxesFontName', fn, ...
           'DefaultLineLineWidth', lw, 'DefaultTextFontSize', fs, ...
           'DefaultTextFontName', fn);

% --- Time parameters ---
dt        = 0.1;
num_steps = 200;
t_sim     = (0:num_steps-1) * dt;

% --- System dynamics (constant-acceleration model) ---
F = [1, dt, 0.5*dt^2; ...
     0,  1,       dt; ...
     0,  0,        1];

f = @(x) F * x;                       % UKF state transition model
h = @(x) [x(1); x(2)];                % Linear measurement function: pos & vel

% --- UKF tuning parameters ---
alpha = 1e-3;
beta  = 2;
kappa = 0;

% --- Initial conditions and noise ---
x_0 = [0; 0; 1];
P   = diag([.1, .1, .1]);          % Initial covariance
Q   = diag([0.001^2, 0.001^2, 0.01^2]);  % Process noise
R   = diag([1^2, 0.1^2]);             % Measurement noise

%% --- Simulation 1: Only Random Process Noise ---
% Create the UKF object
ukf_noise = UKF(f, h, x_0, P, Q, R, alpha, beta, kappa);

% Generate true state and noisy measurements for noise-only case
true_state_noise = zeros(3, num_steps);
meas            = zeros(2, num_steps);
true_state_noise(:, 1) = x_0;

for k = 1:(num_steps-1)
    % True state evolution with only random process noise
    true_state_noise(:, k+1) = f(true_state_noise(:, k)) + mvnrnd(zeros(3,1), Q)';
    % Generate noisy measurement using the linear function
    meas(:, k+1) = h(true_state_noise(:, k+1)) + mvnrnd([0;0], R)';
end

% Run the UKF for the noise-only case
est_state_noise = zeros(3, num_steps);
cov_noise       = zeros(3, 3, num_steps); % Store covariance
est_state_noise(:, 1) = x_0;
cov_noise(:, :, 1) = P;

for k = 2:num_steps
    ukf_noise = predict(ukf_noise);
    ukf_noise = correct(ukf_noise, meas(:, k));
    est_state_noise(:, k) = ukf_noise.getState();
    cov_noise(:, :, k)    = ukf_noise.getCovariance();
end

%% --- Simulation 2: Random Noise + Sinusoidal Perturbation ---
% Reset and create a new UKF object for the second run
ukf_pert = UKF(f, h, x_0, P, Q, R, alpha, beta, kappa);

% Generate true state and noisy measurements for the perturbation case
true_state_pert = zeros(3, num_steps);
meas_pert       = zeros(2, num_steps);
true_state_pert(:, 1) = x_0;

% Define the frequency and acceleration amplitude
freq    = 2*pi/10;
acc_amp = freq;

% Generate pure perturbation states for plotting (analytic solution)
pert_state   = zeros(3, num_steps);
pert_state(:, 1) = x_0;
vel_amplitude = acc_amp / freq;
pos_amplitude = vel_amplitude / freq;

for k = 1:(num_steps-1)
    % The deterministic sinusoidal acceleration input
    t = t_sim(k);
    perturbation_accel = acc_amp * cos(freq * t);

    % Store the pure perturbation states separately for plotting
    t_next    = t_sim(k+1);
    pert_vel  = vel_amplitude * sin(freq * t_next);
    pert_pos  = -pos_amplitude * cos(freq * t_next) + pos_amplitude; % start @ 0
    pert_acc_plot = acc_amp * cos(freq * t_next);
    pert_state(:, k+1) = [pert_pos; pert_vel; pert_acc_plot];

    % Correctly update the true state with the sinusoidal acceleration
    current_true_x = true_state_pert(:, k);
    new_state = [current_true_x(1) + current_true_x(2)*dt + 0.5*perturbation_accel*dt^2;
                 current_true_x(2) + perturbation_accel*dt;
                 perturbation_accel];
    % Update the true state with sinusoidal acceleration and random process noise
    true_state_pert(:, k+1) = new_state + mvnrnd(zeros(3,1), Q)';

    % Generate noisy measurement using the linear function
    meas_pert(:, k+1) = h(true_state_pert(:, k+1)) + mvnrnd([0;0], R)';
end

% Run the UKF for the perturbation case
est_state_pert = zeros(3, num_steps);
cov_pert       = zeros(3, 3, num_steps); % Store covariance
est_state_pert(:, 1) = x_0;
cov_pert(:, :, 1) = P;

for k = 2:num_steps
    ukf_pert = predict(ukf_pert);
    ukf_pert = correct(ukf_pert, meas_pert(:, k));
    est_state_pert(:, k) = ukf_pert.getState();
    cov_pert(:, :, k)    = ukf_pert.getCovariance();
end

%% --- Simulation 3: UKF with DMC for Sinusoidal Perturbation Case ---
% Reset initial conditions for the DMC case
% The augmented state is [pos; vel; accel; dmc_accel_state]
% dmc_accel_state is the Singer model state representing the deviation
x_0_dmc = [0; 0; 0; 0];
P_dmc   = diag([.1, .1, .1, .1]);   % 4x4 covariance for the augmented state

% DMC parameters
dmc_beta     = 0.005;
dmc_sigma_u2 = 0.2456;

% Define the state transition function for the augmented DMC state.
f_dmc = create_f_dmc(dmc_beta);

function f_dmc_handle = create_f_dmc(dmc_beta)
    % Discrete-time closed-form propagation for CA state under an exponentially
    % decaying excitation a_d (Singer): a_d(k+1) = exp(-beta*dt) a_d(k)
    % x = [p; v; a; a_d]
    f_dmc_handle = @(x_aug, dt) [ ...
        x_aug(1) + x_aug(2)*dt + x_aug(3)*0.5*dt^2 ...
            + x_aug(4)/dmc_beta^2 * (dmc_beta*dt - 1 + exp(-dmc_beta*dt));    % Position
        x_aug(2) + x_aug(3)*dt ...
            + x_aug(4)/dmc_beta * (1 - exp(-dmc_beta*dt));                     % Velocity
        x_aug(3) + x_aug(4)/dmc_beta * (1 - exp(-dmc_beta*dt));                % Acceleration
        x_aug(4) * exp(-dmc_beta*dt)];                                         % DMC accel
end

% Measurement function for the augmented state (measure pos and vel)
h_dmc = @(x_aug) [x_aug(1); x_aug(2)];

% Base process noise for the non-DMC states (assume small noise here)
Q_base = diag([0.001^2, 0.001^2, 0.001^2]);

% Create the UKF_DMC object
ukf_dmc = UKF_DMC(f_dmc, h_dmc, x_0_dmc, P_dmc, dmc_beta, dmc_sigma_u2, ...
                  Q_base, R, alpha, beta, kappa, t_sim(1));

% Run the UKF with DMC
est_state_dmc = zeros(4, num_steps);
cov_dmc       = zeros(4, 4, num_steps);
est_state_dmc(:, 1) = x_0_dmc;
cov_dmc(:, :, 1)    = P_dmc;

for k = 2:num_steps
    ukf_dmc = predict(ukf_dmc, t_sim(k));
    ukf_dmc = correct(ukf_dmc, meas_pert(:, k));
    est_state_dmc(:, k) = ukf_dmc.getState();
    cov_dmc(:, :, k)    = ukf_dmc.getCovariance();
end

%% --- Plotting Main Figure (3x3 Subplots with single legend) ---
t_plot = 1:num_steps;
figure;

% --- Column 1: Noise only ---
subplot(3, 3, 1);
h1 = plot(t_plot, true_state_noise(1, :), 'k-'); hold on;
h2 = plot(t_plot, meas(1, :), 'ro', 'MarkerSize', 3);
h3 = plot(t_plot, est_state_noise(1, :), 'b--');
title('UKF');
ylabel('Position'); grid on;
legend([h1, h2, h3], 'True', 'Measurements', 'Estimate', 'Location', 'northwest');

subplot(3, 3, 4);
plot(t_plot, true_state_noise(2, :), 'k-'); hold on;
plot(t_plot, meas(2, :), 'ro', 'MarkerSize', 3);
plot(t_plot, est_state_noise(2, :), 'b--');
ylabel('Velocity'); grid on;

subplot(3, 3, 7);
plot(t_plot, true_state_noise(3, :), 'k-'); hold on;
plot(t_plot, est_state_noise(3, :), 'b--');
xlabel('Time Step'); ylabel('Acceleration'); grid on;

% --- Column 2: Noise + Sinusoidal Perturbation (Standard UKF) ---
subplot(3, 3, 2);
plot(t_plot, true_state_pert(1, :), 'k-'); hold on;
plot(t_plot, meas_pert(1, :), 'ro', 'MarkerSize', 3);
plot(t_plot, est_state_pert(1, :), 'b--');
title('UKF with Perturbation'); grid on;

subplot(3, 3, 5);
plot(t_plot, true_state_pert(2, :), 'k-'); hold on;
plot(t_plot, meas_pert(2, :), 'ro', 'MarkerSize', 3);
plot(t_plot, est_state_pert(2, :), 'b--'); grid on;

subplot(3, 3, 8);
plot(t_plot, true_state_pert(3, :), 'k-'); hold on;
plot(t_plot, est_state_pert(3, :), 'b--');
xlabel('Time Step'); grid on;

% --- Column 3: Noise + Sinusoidal Perturbation (DMC-UKF) ---
subplot(3, 3, 3);
plot(t_plot, true_state_pert(1, :), 'k-'); hold on;
plot(t_plot, meas_pert(1, :), 'ro', 'MarkerSize', 3);
plot(t_plot, est_state_dmc(1, :) + est_state_dmc(4, :), 'b--');  % pos + DMC contribution
title('UKF-DMC with Perturbation'); grid on;

subplot(3, 3, 6);
plot(t_plot, true_state_pert(2, :), 'k-'); hold on;
plot(t_plot, meas_pert(2, :), 'ro', 'MarkerSize', 3);
plot(t_plot, est_state_dmc(2, :), 'b--'); grid on;

subplot(3, 3, 9);
plot(t_plot, true_state_pert(3, :), 'k-'); hold on;
plot(t_plot, est_state_dmc(3, :) + est_state_dmc(4, :), 'b--');
xlabel('Time Step'); grid on;

%% --- Plotting Errors Figure (3x3 Subplots with 3-sigma bounds) ---
% Calculate State Errors
error_noise      = true_state_noise - est_state_noise;
error_with_pert  = true_state_pert  - est_state_pert;

% For DMC, compare true [pos; vel; acc] against estimated [pos; vel; acc+acc_dmc]
est_pos_dmc_tot   = est_state_dmc(1, :);
est_vel_dmc_tot   = est_state_dmc(2, :);
est_acc_dmc_tot   = est_state_dmc(3, :) + est_state_dmc(4, :);
est_state_dmc_tot = [est_pos_dmc_tot; est_vel_dmc_tot; est_acc_dmc_tot];

error_dmc = true_state_pert - est_state_dmc_tot;

figure;

% Get standard deviations (±3σ bands) from covariance matrices
sig_noise_only = zeros(3, num_steps);
sig_pert       = zeros(3, num_steps);
sigma_dmc      = zeros(3, num_steps);

% Transform augmented covariance P_aug (4x4) -> total 3-state covariance
% [pos; vel; acc_tot], where acc_tot = acc + a_dmc
T = [1 0 0 0;   % pos
     0 1 0 0;   % vel
     0 0 1 1];  % acc_tot

for k = 1:num_steps
    % Standard UKF sigmas
    sig_noise_only(:, k) = sqrt(diag(cov_noise(:, :, k)));
    sig_pert(:, k)       = sqrt(diag(cov_pert(:, :, k)));

    % DMC-UKF: project augmented covariance to total-state space
    Ck    = cov_dmc(:, :, k);      % 4x4 covariance of [pos; vel; acc; a_dmc]
    P_tot = T * Ck * T.';          % 3x3 covariance of [pos; vel; acc_tot]
    sigma_dmc(:, k) = sqrt(diag(P_tot));
end

% Helper for shaded ±3σ patch
xrow = t_plot(:).';  % ensure row vector for patch x-coords

% Column 1: Errors for Noise Only
subplot(3, 3, 1);
hold on;
patch([xrow, fliplr(xrow)], [3*sig_noise_only(1, :), fliplr(-3*sig_noise_only(1, :))], ...
      'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(t_plot, error_noise(1, :), 'r-');
title('Position Error (UKF)'); grid on;

subplot(3, 3, 4);
hold on;
patch([xrow, fliplr(xrow)], [3*sig_noise_only(2, :), fliplr(-3*sig_noise_only(2, :))], ...
      'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(t_plot, error_noise(2, :), 'r-');
title('Velocity Error (UKF)'); grid on;

subplot(3, 3, 7);
hold on;
patch([xrow, fliplr(xrow)], [3*sig_noise_only(3, :), fliplr(-3*sig_noise_only(3, :))], ...
      'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(t_plot, error_noise(3, :), 'r-');
xlabel('Time Step'); title('Acceleration Error (UKF)'); grid on;

% Column 2: Errors for Noise + Sinusoidal Perturbation (Standard UKF)
subplot(3, 3, 2);
hold on;
patch([xrow, fliplr(xrow)], [3*sig_pert(1, :), fliplr(-3*sig_pert(1, :))], ...
      'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(t_plot, error_with_pert(1, :), 'r-');
title('Position Error (UKF)'); grid on;

subplot(3, 3, 5);
hold on;
patch([xrow, fliplr(xrow)], [3*sig_pert(2, :), fliplr(-3*sig_pert(2, :))], ...
      'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(t_plot, error_with_pert(2, :), 'r-');
title('Velocity Error (UKF)'); grid on;

subplot(3, 3, 8);
hold on;
patch([xrow, fliplr(xrow)], [3*sig_pert(3, :), fliplr(-3*sig_pert(3, :))], ...
      'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(t_plot, error_with_pert(3, :), 'r-');
xlabel('Time Step'); title('Acceleration Error (UKF)'); grid on;

% Column 3: Errors for Noise + Sinusoidal Perturbation (DMC-UKF)
subplot(3, 3, 3);
hold on;
patch([xrow, fliplr(xrow)], [3*sigma_dmc(1, :), fliplr(-3*sigma_dmc(1, :))], ...
      'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(t_plot, error_dmc(1, :), 'r-');
title('Position Error (DMC-UKF)'); grid on;

subplot(3, 3, 6);
hold on;
patch([xrow, fliplr(xrow)], [3*sigma_dmc(2, :), fliplr(-3*sigma_dmc(2, :))], ...
      'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(t_plot, error_dmc(2, :), 'r-');
title('Velocity Error (DMC-UKF)'); grid on;

subplot(3, 3, 9);
hold on;
patch([xrow, fliplr(xrow)], [3*sigma_dmc(3, :), fliplr(-3*sigma_dmc(3, :))], ...
      'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(t_plot, error_dmc(3, :), 'r-');
xlabel('Time Step'); title('Acceleration Error (DMC-UKF)'); grid on;

% %% --- Plotting the pure perturbation components ---
% figure;
% t_perturb = (0:num_steps-1) * dt;
% 
% subplot(3, 3, 1);
% plot(t_perturb, pert_state(1, :), 'k-');
% title('Sinusoidal Position Perturbation');
% ylabel('Position (m)'); xlabel('Time (s)'); grid on;
% 
% subplot(3, 3, 2);
% plot(t_perturb, pert_state(2, :), 'k-');
% title('Sinusoidal Velocity Perturbation');
% ylabel('Velocity (m/s)'); xlabel('Time (s)'); grid on;
% 
% subplot(3, 3, 3);
% plot(t_perturb, pert_state(3, :), 'k-'); hold on;
% % Plot the DMC-estimated acceleration perturbation (4th row)
% plot(t_plot * dt, est_state_dmc(4, :), 'b--');
% title('Sinusoidal Acceleration Perturbation');
% ylabel('Acceleration (m/s^2)'); xlabel('Time (s)');
% legend('True Perturbation', 'UKF-DMC Estimate', 'Location', 'best');
% grid on;
% 
% sgtitle('Individual Sinusoidal Perturbation Components Over Time');