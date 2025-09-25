%% spacecraft_with_UKF_2figs_cols.m
clear; clc; close all;

%% Plot style
style = PlotStyle(); style.apply();   % fs=12, fn='Times New Roman', lw=2

%% Spacecraft dynamics
J = diag([10 8 5]);

% Open-loop torque (before controller turns on): steady + short impact pulse
tau_open = @(t,x) [0.05; -0.02; 0.03] + (t>=0.1 && t<=0.6)*[2; -1; 1];

% Build spacecraft with a placeholder torqueFcn; we will overwrite each step
sc = Spacecraft(J, tau_open);

% Desired target attitude (scalar-first quaternion)
qd = [1;0;0;0];

% PD controller (use AttitudePD class)
ctl = AttitudePD('AutoTune',true,'J',J,'OmegaN',0.8,'Zeta',0.9,'TauMax',[2;2;2]);
ctl = ctl.setDesired(qd);

x0_true = [0;0;0; 1;0;0;0];    % [w; q]

%% UKF setup
dt    = 0.1;
t_end = 20;
t     = 0:dt:t_end;
N     = numel(t);

% One-step transition for UKF (internal RK4 over fixed dt)
fFcn  = @(x,u) rk4(@(x) sc.attitudeDynamics(0,x), x, dt);
hFcn  = @(x) x;

x0_est = x0_true + [0.05; -0.03; 0.04; 0; 0; 0; 0];

Q = diag([5e-5*ones(1,3), 1e-6*ones(1,4)]);
R = diag([(0.01)^2*ones(1,3), (0.01)^2*ones(1,4)]);

ukf = UKF(fFcn, hFcn, x0_est, 1e-2*eye(7), Q, R, 1e-3, 2, 0);

%% Simulation loop (truth, noisy measurements, UKF, controller switching)
Xtrue = zeros(7,N); Xtrue(:,1) = x0_true;
Xest  = zeros(7,N); Xest(:,1)  = x0_est;
Zmeas = zeros(7,N); Zmeas(:,1) = x0_true;
P_hist = zeros(7,7,N); P_hist(:,:,1) = ukf.getCovariance();

% Measurement noise
sigma_w      = 0.03;          % rad/s noise
sigma_q_rad  = 1.5*pi/180;    % quaternion multiplicative noise (~1.5 deg)

% Controller switching logic (estimate-based, NO TIMER)
theta_on_deg   = 5;                  % attitude error angle threshold (deg)
theta_on       = theta_on_deg*pi/180;
w_on           = 0.05;               % rate magnitude threshold (rad/s)
controlOn      = false;              % latched once true
onCount        = 0;                  % debounce counter
onCountNeeded  = 2;                  % need this many consecutive steps

for k = 2:N
    % ====== 1) Decide control torque for next step (dynamics-only trigger) ======
    xhat_prev = Xest(:,k-1);
    qhat = xhat_prev(4:7); qhat = qhat/norm(qhat);
    what = xhat_prev(1:3);

    % Quaternion error q_e = qd^{-1} ⊗ qhat (scalar-first), shortest rotation
    qd0 = qd(1); qdv = qd(2:4);
    qinv = [qd0; -qdv];
    s1=qinv(1); v1=qinv(2:4);
    s2=qhat(1); v2=qhat(2:4);
    q_e0 = s1*s2 - dot(v1,v2);
    q_ev = s1*v2 + s2*v1 + cross(v1,v2);
    if q_e0 < 0, q_e0 = -q_e0; q_ev = -q_ev; end
    theta_err = 2*acos( max(-1,min(1,q_e0)) );   % attitude error angle (rad)

    % Debounced dynamics-based trigger
    if (theta_err > theta_on) || (norm(what) > w_on)
        onCount = onCount + 1;
    else
        onCount = 0;
    end
    if ~controlOn && (onCount >= onCountNeeded)
        controlOn = true;   % latch ON
    end

    % Select torque for this step
    if controlOn
        % Closed-loop: PD torque from estimated state (qhat, what)
        tau_cmd = ctl.torque(qhat, what);
        sc.torqueFcn = @(tt,xx) tau_cmd;  % hold constant over this dt
    else
        % Open-loop torque (with impact pulse)
        sc.torqueFcn = tau_open;
    end

    % ====== 2) Propagate TRUE state over [t(k-1), t(k)] with chosen torque ======
    xk = rk4(@(x) sc.attitudeDynamics(t(k-1), x), Xtrue(:,k-1), dt);
    xk(4:7) = xk(4:7)/norm(xk(4:7));
    Xtrue(:,k) = xk;

    % ====== 3) Generate noisy measurement z_k = [w_meas; q_meas] ======
    % Angular velocity: additive Gaussian
    w_meas = xk(1:3) + sigma_w*randn(3,1);

    % Quaternion: multiplicative small-angle noise (inlined)
    eta = sigma_q_rad * randn(3,1); th = norm(eta);
    if th < 1e-12
        dq = [1;0;0;0];
    else
        ax = eta/th; half = 0.5*th;
        dq = [cos(half); ax*sin(half)];
    end
    % q_meas = dq ⊗ q_true (scalar-first multiply, inline)
    s1 = dq(1); v1 = dq(2:4);
    s2 = xk(4); v2 = xk(5:7);
    s  = s1*s2 - dot(v1,v2);
    v  = s1*v2 + s2*v1 + cross(v1,v2);
    q_meas = [s; v]; q_meas = q_meas / norm(q_meas);

    Zmeas(:,k) = [w_meas; q_meas];

    % ====== 4) UKF predict/correct ======
    ukf = ukf.predict();
    ukf = ukf.correct(Zmeas(:,k));
    xhat = ukf.getState(); xhat(4:7) = xhat(4:7)/norm(xhat(4:7));
    Xest(:,k) = xhat;
    P_hist(:,:,k) = ukf.getCovariance();
end

%% FIGURE 1: Quaternion (4 rows × 2 cols: State | Error ±3σ)
figure('Name','Quaternion: States and Errors with 3\sigma');
names_q = {'q_0','q_1','q_2','q_3'};
for i = 1:4
    idx = 3 + i;                 % 4..7
    err = Xtrue(idx,:) - Xest(idx,:);
    sig = sqrt(squeeze(P_hist(idx,idx,:))).';

    % Left column: state
    subplot(4,2,2*i-1);
    plot(t, Xtrue(idx,:), 'k-'); hold on;
    plot(t, Zmeas(idx,:), 'ro', 'MarkerSize', 3);
    plot(t, Xest(idx,:),  'b--');
    if i==1, title('State'); end
    ylabel(names_q{i}); grid on;
    if i==4, xlabel('Time (s)'); end
    if i==1, legend('true','meas','ukf','Location','best'); end

    % Right column: error ±3σ centered at 0
    subplot(4,2,2*i);
    patch([t, fliplr(t)], [ 3*sig, fliplr(-3*sig) ], ...
          'b', 'FaceAlpha', 0.2, 'EdgeColor','none'); hold on;
    plot(t, err, 'r-');
    yline(0,'k-');
    if i==1, title('Error & 3\sigma'); end
    grid on;
    if i==4, xlabel('Time (s)'); end
end

%% FIGURE 2: Body rates (3 rows × 2 cols: State | Error ±3σ)
figure('Name','Body Rates: States and Errors with 3\sigma');
names_w = {'\omega_x','\omega_y','\omega_z'};
for i = 1:3
    idx = i;                      % 1..3
    err = Xtrue(idx,:) - Xest(idx,:);
    sig = sqrt(squeeze(P_hist(idx,idx,:))).';

    % Left column: state
    subplot(3,2,2*i-1);
    plot(t, Xtrue(idx,:), 'k-'); hold on;
    plot(t, Zmeas(idx,:), 'ro', 'MarkerSize', 3);
    plot(t, Xest(idx,:),  'b--');
    if i==1, title('State'); end
    ylabel(names_w{i}); grid on;
    if i==3, xlabel('Time (s)'); end
    if i==1, legend('true','meas','ukf','Location','best'); end

    % Right column: error ±3σ
    subplot(3,2,2*i);
    patch([t, fliplr(t)], [ 3*sig, fliplr(-3*sig) ], ...
          'b', 'FaceAlpha', 0.2, 'EdgeColor','none'); hold on;
    plot(t, err, 'r-');
    yline(0,'k-');
    if i==1, title('Error & 3\sigma'); end
    grid on;
    if i==3, xlabel('Time (s)'); end
end

%% Local integrator
function x_next = rk4(f, x, h)
    k1=f(x); k2=f(x+0.5*h*k1); k3=f(x+0.5*h*k2); k4=f(x+h*k3);
    x_next = x + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
end