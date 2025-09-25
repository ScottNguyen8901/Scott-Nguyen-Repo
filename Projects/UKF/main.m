%% dmc_ukf.m  —  Demo/driver for UKF vs. UKF_DMC on a constant-accel model
% Scott-style: tidy headers, reproducible RNG, and Times New Roman plots.

%% Setup
close all; clear; clc;
rng(1);  % reproducibility

% --- Default plotting properties ---
fs = 12;             % default font size
fn = 'Times New Roman';
lw = 2;              % line width
set(groot, 'DefaultAxesFontSize', fs, ...
           'DefaultAxesFontName', fn, ...
           'DefaultLineLineWidth', lw, ...
           'DefaultTextFontSize', fs, ...
           'DefaultTextFontName', fn);

%% Time grid
dt       = 0.1;          % s  (10 Hz)
num_steps = 200;
t_sim    = (0:num_steps-1) * dt;

%% System (original) state: [position; velocity; acceleration]
% Discrete-time kinematics for constant-acceleration (used for the UKF model)
f_disc = @(x,dt) [ x(1) + x(2)*dt + 0.5*x(3)*dt^2; ...
                   x(2) + x(3)*dt; ...
                   x(3) ];
% Measurement: we observe position and velocity
h_meas = @(x) x([1,2]);

%% UKF tuning
alpha = 1e-3;   beta = 2;   kappa = 0;

%% Initial conditions & noise
initial_x = [0; 0; 0];                    % x0
P0        = diag([0.1, 0.1, 0.1]);        % P0

Q         = diag([1e-3^2, 1e-3^2, 1e-2^2]);  % process noise (model)
R         = diag([1^2, 0.1^2]);              % measurement noise (pos, vel)

%% ---------- Simulation 1: Noise only ------------------------------------
% Create the standard UKF (uses dt as the "u" argument)
f_std = @(x,u) f_disc(x, u);         % pass dt via u
ukf_noise_only = UKF(f_std, h_meas, initial_x, P0, Q, R, alpha, beta, kappa);

true_state_noise_only = zeros(3, num_steps);
meas_noise_only       = zeros(2, num_steps);

true_state_noise_only(:,1) = initial_x;
meas_noise_only(:,1)       = h_meas(true_state_noise_only(:,1)) + mvnrnd([0,0], R).';

for k = 1:num_steps-1
    % propagate true with process noise
    xk1 = f_disc(true_state_noise_only(:,k), dt) + mvnrnd([0 0 0], Q).';
    true_state_noise_only(:,k+1) = xk1;

    % measurement
    meas_noise_only(:,k+1) = h_meas(xk1) + mvnrnd([0 0], R).';
end

% Run the UKF (noise-only case)
est_noise_only = zeros(3, num_steps);
cov_noise_only = zeros(3, 3, num_steps);
est_noise_only(:,1) = initial_x;   cov_noise_only(:,:,1) = P0;

for k = 2:num_steps
    ukf_noise_only = ukf_noise_only.predict(dt);
    ukf_noise_only = ukf_noise_only.correct(meas_noise_only(:,k));
    est_noise_only(:,k) = ukf_noise_only.getState();
    cov_noise_only(:,:,k) = ukf_noise_only.getCovariance();
end

%% ---------- Simulation 2: Noise + sinusoidal perturbation (standard UKF) -
% Recreate a fresh standard UKF
ukf_std = UKF(f_std, h_meas, initial_x, P0, Q, R, alpha, beta, kappa);

true_state_w = zeros(3, num_steps);
meas_w       = zeros(2, num_steps);
true_state_w(:,1) = initial_x;
meas_w(:,1)  = h_meas(true_state_w(:,1)) + mvnrnd([0 0], R).';

% Perturbation a_p(t) = A*cos(2π f t)
freq = 2*pi/10;       % rad/s  (period ~10 steps)
Aacc = freq;          % choose amplitude ~ freq for nice dynamics

% Analytic “pure perturbation” components for plotting (pos, vel, accel)
perturb_state = zeros(3, num_steps);  % [p; v; a] caused by sinusoid alone
pos_amp  = Aacc / freq^2;
vel_amp  = Aacc / freq;
for k = 1:num_steps-1
    t = t_sim(k); t_next = t_sim(k+1);
    a_p   = Aacc * cos(freq*t);
    a_p_n = Aacc * cos(freq*t_next);

    % store pure perturbation components
    perturb_state(:,k) = [pos_amp*(1-cos(freq*t));
                          vel_amp*sin(freq*t);
                          a_p];

    % update true state with sinusoidal acceleration + process noise
    xk     = true_state_w(:,k);
    xk1    = [xk(1) + xk(2)*dt + 0.5*(xk(3) + a_p)*dt^2;
              xk(2) + (xk(3) + a_p)*dt;
              xk(3)] + mvnrnd([0 0 0], Q).';
    true_state_w(:,k+1) = xk1;

    meas_w(:,k+1) = h_meas(xk1) + mvnrnd([0 0], R).';
end
perturb_state(:,end) = [pos_amp*(1-cos(freq*t_sim(end)));
                        vel_amp*sin(freq*t_sim(end));
                        Aacc*cos(freq*t_sim(end))];

% Run the UKF on the perturbed data
est_w = zeros(3, num_steps);
cov_w = zeros(3,3,num_steps);
est_w(:,1) = initial_x; cov_w(:,:,1) = P0;

for k = 2:num_steps
    ukf_std = ukf_std.predict(dt);
    ukf_std = ukf_std.correct(meas_w(:,k));
    est_w(:,k) = ukf_std.getState();
    cov_w(:,:,k) = ukf_std.getCovariance();
end

%% ---------- Simulation 3: DMC-UKF for the same perturbation -------------
% Augmented state: [pos; vel; acc; dmc_pos; dmc_vel; dmc_acc_state]
initial_x_dmc = [initial_x; 0; 0; 0];
P0_dmc        = blkdiag(P0, diag([0.1, 0.1, 0.1]));

% Base transition for original 3 states
f_base = @(x,dt) [ x(1) + x(2)*dt + 0.5*x(3)*dt^2;
                   x(2) + x(3)*dt;
                   x(3) ];

% DMC parameters
dmc_beta    = freq;    % GM rate ~= sinusoid freq works well
dmc_sigma_u2 = 1e-4;   % driving noise variance

% Build the augmented transition using a nested creator (captures beta)
f_dmc = create_f_dmc(dmc_beta);

% Measurement of augmented: we still measure position & velocity (base only)
h_dmc = @(x_aug) x_aug([1,2]);

% Process noise for non-DMC part (base states) — small; DMC Q is handled inside UKF_DMC
Q_base = diag([1e-3^2, 1e-3^2, 1e-2^2]);

% Create the DMC-UKF
ukf_dmc = UKF_DMC(f_dmc, h_dmc, initial_x_dmc, P0_dmc, ...
                  Q_base, R, alpha, beta, kappa, t_sim(1), dmc_beta, dmc_sigma_u2);

% Run the DMC-UKF on the same perturbed measurements
est_dmc = zeros(6, num_steps);
cov_dmc = zeros(6, 6, num_steps);
est_dmc(:,1)   = initial_x_dmc;
cov_dmc(:,:,1) = P0_dmc;

for k = 2:num_steps
    ukf_dmc = ukf_dmc.predict(t_sim(k));                 % pass absolute time
    ukf_dmc = ukf_dmc.correct(meas_w(:,k));
    est_dmc(:,k)   = ukf_dmc.getState();
    cov_dmc(:,:,k) = ukf_dmc.getCovariance();
end

%% ---------- PLOTTING: States (3×3) --------------------------------------
t_plot = 1:num_steps;

figure;
% Column 1 — Noise only
subplot(3,3,1);
plot(t_plot, true_state_noise_only(1,:), 'k-'); hold on;
plot(t_plot, meas_noise_only(1,:), 'ro', 'MarkerSize', 3);
plot(t_plot, est_noise_only(1,:), 'b--');
title('Position'); ylabel('Position'); grid on;

subplot(3,3,4);
plot(t_plot, true_state_noise_only(2,:), 'k-'); hold on;
plot(t_plot, meas_noise_only(2,:), 'ro', 'MarkerSize', 3);
plot(t_plot, est_noise_only(2,:), 'b--');
title('Velocity'); ylabel('Velocity'); grid on;

subplot(3,3,7);
plot(t_plot, true_state_noise_only(3,:), 'k-'); hold on;
plot(t_plot, est_noise_only(3,:), 'b--');
title('Acceleration'); xlabel('Time Step'); ylabel('Acceleration'); grid on;

% Column 2 — Noise + sinusoid (standard UKF)
subplot(3,3,2);
plot(t_plot, true_state_w(1,:), 'k-'); hold on;
plot(t_plot, meas_w(1,:), 'ro', 'MarkerSize', 3);
plot(t_plot, est_w(1,:), 'b--');
title('Position (UKF)'); grid on;

subplot(3,3,5);
plot(t_plot, true_state_w(2,:), 'k-'); hold on;
plot(t_plot, meas_w(2,:), 'ro', 'MarkerSize', 3);
plot(t_plot, est_w(2,:), 'b--');
title('Velocity (UKF)'); grid on;

subplot(3,3,8);
plot(t_plot, true_state_w(3,:), 'k-'); hold on;
plot(t_plot, est_w(3,:), 'b--');
title('Acceleration (UKF)'); xlabel('Time Step'); grid on;

% Column 3 — Noise + sinusoid (DMC-UKF). Plot TOTAL = base + DMC for p,v,a
subplot(3,3,3);
plot(t_plot, true_state_w(1,:), 'k-'); hold on;
plot(t_plot, meas_w(1,:), 'ro', 'MarkerSize', 3);
plot(t_plot, est_dmc(1,:)+est_dmc(4,:), 'b--');
title('Position (DMC-UKF)'); grid on;

subplot(3,3,6);
plot(t_plot, true_state_w(2,:), 'k-'); hold on;
plot(t_plot, meas_w(2,:), 'ro', 'MarkerSize', 3);
plot(t_plot, est_dmc(2,:)+est_dmc(5,:), 'b--');
title('Velocity (DMC-UKF)'); grid on;

subplot(3,3,9);
plot(t_plot, true_state_w(3,:), 'k-'); hold on;
plot(t_plot, est_dmc(3,:)+est_dmc(6,:), 'b--');
title('Acceleration (DMC-UKF)'); xlabel('Time Step'); grid on;

annotation('textbox',[0.13, 0.95, 0.2, 0.04],'String','Noise Only', ...
    'EdgeColor','none','HorizontalAlignment','center','FontWeight','bold');
annotation('textbox',[0.41, 0.95, 0.2, 0.04],'String','Noise + Perturbation (UKF)', ...
    'EdgeColor','none','HorizontalAlignment','center','FontWeight','bold');
annotation('textbox',[0.69, 0.95, 0.2, 0.04],'String','Noise + Perturbation (DMC-UKF)', ...
    'EdgeColor','none','HorizontalAlignment','center','FontWeight','bold');

%% ---------- PLOTTING: Errors with 3-sigma bounds (3×3) ------------------
err_noise   = true_state_noise_only - est_noise_only;
err_w       = true_state_w - est_w;

% standard deviations
sigma_noise = zeros(3, num_steps);
sigma_w     = zeros(3, num_steps);
sigma_dmc   = zeros(3, num_steps);  % combine base+DMC covariance

for k = 1:num_steps
    sigma_noise(:,k) = sqrt(diag(cov_noise_only(:,:,k)));
    sigma_w(:,k)     = sqrt(diag(cov_w(:,:,k)));

    % Combine base and DMC covariance blocks for total p/v/a
    C = cov_dmc(:,:,k);
    % totals: [1+4, 2+5, 3+6]
    C_total = [ C(1,1)+C(4,4)+2*C(1,4),   C(1,2)+C(4,5)+C(1,5)+C(4,2),   C(1,3)+C(4,6)+C(1,6)+C(4,3);  % p with (v,a)
                C(2,1)+C(5,4)+C(2,4)+C(5,1), C(2,2)+C(5,5)+2*C(2,5),     C(2,3)+C(5,6)+C(2,6)+C(5,3);
                C(3,1)+C(6,4)+C(3,4)+C(6,1), C(3,2)+C(6,5)+C(3,5)+C(6,2), C(3,3)+C(6,6)+2*C(3,6) ];
    sigma_dmc(:,k) = sqrt(diag(C_total));
end

err_w_dmc = [ true_state_w(1,:) - (est_dmc(1,:)+est_dmc(4,:)); ...
              true_state_w(2,:) - (est_dmc(2,:)+est_dmc(5,:)); ...
              true_state_w(3,:) - (est_dmc(3,:)+est_dmc(6,:)) ];

figure;

% Column 1 — Noise only errors
for r = 1:3
    subplot(3,3,3*(r-1)+1); hold on;
    patch([t_plot, fliplr(t_plot)], [ 3*sigma_noise(r,:), fliplr(-3*sigma_noise(r,:)) ], ...
          'b', 'FaceAlpha', 0.2, 'EdgeColor','none');
    plot(t_plot, err_noise(r,:), 'r-');
    if r==1, title('Position Error'); elseif r==2, title('Velocity Error'); else, title('Acceleration Error'); end
    grid on;
end

% Column 2 — UKF errors with sinusoid
for r = 1:3
    subplot(3,3,3*(r-1)+2); hold on;
    patch([t_plot, fliplr(t_plot)], [ 3*sigma_w(r,:), fliplr(-3*sigma_w(r,:)) ], ...
          'b', 'FaceAlpha', 0.2, 'EdgeColor','none');
    plot(t_plot, err_w(r,:), 'r-');
    if r==1, title('Position Error (UKF)'); elseif r==2, title('Velocity Error (UKF)'); else, title('Acceleration Error (UKF)'); end
    grid on;
end

% Column 3 — DMC-UKF errors with sinusoid
for r = 1:3
    subplot(3,3,3*(r-1)+3); hold on;
    patch([t_plot, fliplr(t_plot)], [ 3*sigma_dmc(r,:), fliplr(-3*sigma_dmc(r,:)) ], ...
          'b', 'FaceAlpha', 0.2, 'EdgeColor','none');
    plot(t_plot, err_w_dmc(r,:), 'r-');
    if r==1, title('Position Error (DMC-UKF)'); elseif r==2, title('Velocity Error (DMC-UKF)'); else, title('Acceleration Error (DMC-UKF)'); end
    grid on;
end

annotation('textbox',[0.13, 0.95, 0.2, 0.04],'String','Noise Only', ...
    'EdgeColor','none','HorizontalAlignment','center','FontWeight','bold');
annotation('textbox',[0.41, 0.95, 0.2, 0.04],'String','Noise + Perturbation (UKF)', ...
    'EdgeColor','none','HorizontalAlignment','center','FontWeight','bold');
annotation('textbox',[0.69, 0.95, 0.2, 0.04],'String','Noise + Perturbation (DMC-UKF)', ...
    'EdgeColor','none','HorizontalAlignment','center','FontWeight','bold');

%% ---------- PLOTTING: Pure perturbation components ----------------------
figure;
subplot(3,1,1);
plot(t_sim, perturb_state(1,:), 'k-');
title('Sinusoidal Position Perturbation'); ylabel('Position (m)'); xlabel('Time (s)'); grid on;

subplot(3,1,2);
plot(t_sim, perturb_state(2,:), 'k-');
title('Sinusoidal Velocity Perturbation'); ylabel('Velocity (m/s)'); xlabel('Time (s)'); grid on;

subplot(3,1,3); hold on;
plot(t_sim, perturb_state(3,:), 'k-');
plot(t_plot, est_dmc(6,:), 'b--');
legend('True Perturbation','UKF-DMC Estimate','Location','best');
title('Sinusoidal Acceleration Perturbation');
ylabel('Acceleration (m/s^2)'); xlabel('Time (s)'); grid on;

sgtitle('Individual Sinusoidal Perturbation Components Over Time');

%% -------------------- Helper: build f_dmc (nested) ----------------------
function f_dmc = create_f_dmc(beta)
% Returns f_dmc(x_aug, dt) for augmented state:
% x_aug = [pos; vel; acc; dmc_pos; dmc_vel; dmc_acc_state]
    f_dmc = @(x_aug, dt) f_dmc_impl(x_aug, dt, beta);
end

function x_next = f_dmc_impl(x_aug, dt, beta)
    x_next = zeros(6,1);

    % Base states (pos/vel/acc) — use base acc + DMC acc as effective accel
    a_eff = x_aug(3) + x_aug(6);
    x_next(1) = x_aug(1) + x_aug(2)*dt + 0.5*a_eff*dt^2;   % position
    x_next(2) = x_aug(2) + a_eff*dt;                        % velocity
    x_next(3) = x_aug(3);                                   % base accel (modeled constant)

    % DMC integrated states
    x_next(4) = x_aug(4) + x_aug(5)*dt + 0.5*x_aug(6)*dt^2;                 % dmc pos
    x_next(5) = x_aug(5) + x_aug(6)*dt;                                     % dmc vel

    % Gauss–Markov (1st order) for dmc_acc_state
    % x6_dot = -beta*x6 + u(t); discretized exactly:
    x_next(6) = x_aug(6)*exp(-beta*dt);                                     % mean prop (noise handled via Q_DMC)
end