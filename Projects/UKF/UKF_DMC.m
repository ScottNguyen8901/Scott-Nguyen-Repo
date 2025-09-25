classdef UKF_DMC
%UKF_DMC  Unscented Kalman filter with Dynamic Model Compensation (DMC).
%
%   This extends a standard (additive-noise) UKF by augmenting the state
%   with a 3-state Gauss–Markov process to capture unmodeled dynamics.
%
%   Models (additive noise):
%       x_{k+1} = f(x_k, dt) + w_k,     w_k ~ N(0, Q_aug(dt))
%       z_k     = h(x_k)     + v_k,     v_k ~ N(0, R)
%
%   Constructor:
%       obj = UKF_DMC(fFcn, hFcn, x0, P0, Q_base, R, ...
%                     alpha, beta_ukf, kappa, t_initial, ...
%                     dmc_beta, dmc_sigma_u2)
%
%       fFcn          : state transition, signature f(x, dt)
%       hFcn          : measurement function, signature h(x)
%       x0, P0        : initial augmented state and covariance
%       Q_base        : process noise for the ORIGINAL part of the state
%                       (size n_orig×n_orig; non-DMC block)
%       R             : measurement noise covariance
%       alpha,beta,kappa : UKF tuning
%       t_initial     : initial timestamp (seconds)
%       dmc_beta      : Gauss–Markov rate parameter β (>0)
%       dmc_sigma_u2  : driving-noise variance σ_u^2 for the DMC process
%
%   Main methods:
%       obj = predict(obj, t_now)
%       obj = correct(obj, z)
%       x   = getState(obj)
%       P   = getCovariance(obj)
%
%   Notes:
%       • Augmented state = [ x_orig ; x_dmc ] with n_dmc = 3 here.
%       • Q_aug is block-diag{ Q_base , Q_DMC(dt) }.
%       • Uses Cholesky with small jitter for numerical robustness.
%

    %====================== Stored data ======================%
    properties (Access = private)
        % Models
        stateTransitionFcn    % f(x, dt)
        measurementFcn        % h(x)

        % Noise / helpers
        Q_base                % process noise for original (non-DMC) block
        R                     % measurement noise
        Q_n_fcn               % @(dt) -> Q_DMC(dt) (3×3)

        % Estimates
        x                     % current state estimate (n×1)
        P                     % current covariance (n×n)

        % Dimensions
        n_orig                % original state dimension
        n_dmc                 % augmented DMC states (fixed = 3)
        n                     % total state dimension (n_orig + n_dmc)
        m                     % measurement dimension

        % UKF scaling
        alpha
        beta_ukf
        kappa
        lambda
        gamma
        Wm                    % 1×(2n+1)
        Wc                    % 1×(2n+1)

        % Time & DMC params
        t_prev                % previous timestamp [s]
        dmc_beta              % β
        dmc_sigma_u2          % σ_u^2
    end

    %====================== Public API ======================%
    methods
        function obj = UKF_DMC(stateTransitionFcn, measurementFcn, ...
                               initial_x, initial_P, Q_base, R, ...
                               alpha, beta_ukf, kappa, t_initial, ...
                               dmc_beta, dmc_sigma_u2)
        %UKF_DMC Constructor. See class header for details.

            % Store models/noise
            obj.stateTransitionFcn = stateTransitionFcn;
            obj.measurementFcn     = measurementFcn;
            obj.Q_base             = Q_base;
            obj.R                  = R;

            % Initial estimates
            obj.x = initial_x(:);
            obj.P = initial_P;

            % DMC parameters
            obj.dmc_beta     = dmc_beta;
            obj.dmc_sigma_u2 = dmc_sigma_u2;

            % Dimensions
            obj.n_orig = numel(initial_x) - 3;   % assuming last 3 are DMC
            obj.n_dmc  = 3;
            obj.n      = obj.n_orig + obj.n_dmc;
            obj.m      = size(R,1);

            % UKF tuning
            obj.alpha   = alpha;
            obj.beta_ukf = beta_ukf;
            obj.kappa   = kappa;

            % Precompute weights / scaling
            obj = obj.calculateWeightsAndScaling();

            % Time
            obj.t_prev = t_initial;

            % DMC process noise handle Q_DMC(dt)
            obj.Q_n_fcn = @(dt) computeDmcProcessNoise(dt, obj.dmc_beta, obj.dmc_sigma_u2);
        end

        function obj = predict(obj, t_now)
        %PREDICT  Propagate sigma points through f(x,dt) and add Q_aug(dt).
            dt = t_now - obj.t_prev;

            % Sigma points from current (x,P)
            Xi = obj.generateSigmaPoints();

            % Propagate through dynamics
            Yi = zeros(obj.n, size(Xi,2));
            for i = 1:size(Xi,2)
                Yi(:,i) = obj.stateTransitionFcn(Xi(:,i), dt);
            end

            % Mean & covariance of prediction
            x_pred = Yi * obj.Wm.';                  % n×1
            P_pred = zeros(obj.n);
            for i = 1:size(Yi,2)
                d = Yi(:,i) - x_pred;
                P_pred = P_pred + obj.Wc(i) * (d * d.');
            end

            % Augmented process noise: block-diag{Q_base, Q_DMC(dt)}
            Q_aug = zeros(obj.n);
            Q_aug(1:obj.n_orig, 1:obj.n_orig) = obj.Q_base;
            Q_aug(obj.n_orig+1:end, obj.n_orig+1:end) = obj.Q_n_fcn(max(dt,0));

            % Commit
            obj.x = x_pred;
            obj.P = obj.symmetrize(P_pred + Q_aug);
            obj.t_prev = t_now;
        end

        function obj = correct(obj, z)
        %CORRECT  Incorporate measurement z via h(x) and R.
            Xi = obj.generateSigmaPoints();          % n×(2n+1)

            % Predicted measurements
            Zi = zeros(obj.m, size(Xi,2));
            for i = 1:size(Xi,2)
                Zi(:,i) = obj.measurementFcn(Xi(:,i));
            end
            z_pred = Zi * obj.Wm.';                  % m×1

            % Innovation covariance S and cross-cov Pxz
            S   = zeros(obj.m);
            Pxz = zeros(obj.n, obj.m);
            for i = 1:size(Xi,2)
                dx = Xi(:,i) - obj.x;
                dz = Zi(:,i) - z_pred;
                S   = S   + obj.Wc(i) * (dz * dz.');
                Pxz = Pxz + obj.Wc(i) * (dx * dz.');
            end
            S = obj.symmetrize(S + obj.R);

            % Gain & update
            K = Pxz / (S + 1e-12*eye(obj.m));
            obj.x = obj.x + K * (z - z_pred);
            obj.P = obj.symmetrize(obj.P - K * S * K.');
        end

        function x_est = getState(obj),      x_est = obj.x; end
        function P_est = getCovariance(obj), P_est = obj.P; end
    end

    %====================== Internals ======================%
    methods (Access = private)
        function obj = calculateWeightsAndScaling(obj)
            obj.lambda = obj.alpha^2 * (obj.n + obj.kappa) - obj.n;
            obj.gamma  = sqrt(obj.n + obj.lambda);

            obj.Wm = zeros(1, 2*obj.n + 1);
            obj.Wc = zeros(1, 2*obj.n + 1);
            obj.Wm(1) = obj.lambda / (obj.n + obj.lambda);
            obj.Wc(1) = obj.Wm(1) + (1 - obj.alpha^2 + obj.beta_ukf);
            obj.Wm(2:end) = 1/(2*(obj.n + obj.lambda));
            obj.Wc(2:end) = obj.Wm(2:end);
        end

        function Xi = generateSigmaPoints(obj)
        %GENERATESIGMAPOINTS  n×(2n+1) matrix of sigma points about (x,P).
            % Robust Cholesky
            jitter = 1e-12;
            ok = false; tries = 0;
            while ~ok && tries < 5
                try
                    A = chol(obj.P + jitter*eye(obj.n), 'lower');
                    ok = true;
                catch
                    jitter = max(jitter*10, 1e-12);
                    tries = tries + 1;
                end
            end
            if ~ok
                [V,D] = eig((obj.P+obj.P')/2);
                D = max(D, 1e-12*eye(obj.n));
                A = chol(V*D*V.', 'lower');
            end

            Xi = zeros(obj.n, 2*obj.n + 1);
            Xi(:,1) = obj.x;
            for i = 1:obj.n
                Xi(:,1+i)        = obj.x + obj.gamma*A(:,i);
                Xi(:,1+obj.n+i)  = obj.x - obj.gamma*A(:,i);
            end
        end

        function M = symmetrize(~, M), M = (M + M.')/2; end
    end
end

%================== DMC process-noise helper ==================%
function Q_n = computeDmcProcessNoise(dt, beta, sigma_u2)
%COMPUTEDMCPROCESSNOISE  Closed-form Q for a 3-state Gauss–Markov chain
% driven by white noise of variance sigma_u2, with rate parameter beta.
%
% This matches the structure visible in your code screenshots: the
% elements below are symmetric with terms in t, exp(-βt), and exp(-2βt).
% (Use your exact continuous-time model if you have one; this is the
% typical closed form used for DMC augmentation.)

t = max(dt, 0);

% Diagonal terms
Q_n_1_1 = sigma_u2 * ( (1/(3*beta^2))*t^3 ...
                     + (1/(beta^3))*t^2 ...
                     + (1/(beta^4))*(1 - 2*exp(-beta*t)) ...
                     + (1/(2*beta^5))*(1 - exp(-2*beta*t)) );

Q_n_2_2 = sigma_u2 * ( (1/(beta^2))*t ...
                     + (2/beta^3)*(1 - exp(-beta*t)) ...
                     + (1/(2*beta^4))*(1 - exp(-2*beta*t)) );

Q_n_3_3 = sigma_u2 * ( (1/(2*beta))*(1 - exp(-2*beta*t)) );

% Off-diagonals (symmetric)
Q_n_1_2 = sigma_u2 * ( (1/(2*beta^2))*t^2 ...
                     + (1/(beta^3))*(1 - exp(-beta*t)) ...
                     + (1/(2*beta^4))*(1 - exp(-2*beta*t)) );

Q_n_1_3 = sigma_u2 * ( (1/(2*beta))*(1 - exp(-2*beta*t)) ...
                     + (1/(beta^2))*exp(-beta*t) );

Q_n_2_3 = sigma_u2 * ( (1/(2*beta^2))*(1 - exp(-2*beta*t)) ...
                     + (1/beta^2)*exp(-beta*t) );

% Assemble symmetric 3×3
Q_n = [Q_n_1_1, Q_n_1_2, Q_n_1_3;
       Q_n_1_2, Q_n_2_2, Q_n_2_3;
       Q_n_1_3, Q_n_2_3, Q_n_3_3];
end