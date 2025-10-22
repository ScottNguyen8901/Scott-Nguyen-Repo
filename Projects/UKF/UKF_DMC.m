classdef UKF_DMC
    %UKF_DMC An unscented Kalman filter class that incorporates Dynamic Model Compensation (DMC).
    %
    % This class extends the standard UKF to handle unmodeled dynamics by
    % augmenting the state vector with a Gauss-Markov process. The state
    % transition function and process noise covariance are specifically
    % designed for this purpose.

    properties (Access = private)
        stateTransitionFcn   % Function handle for the nonlinear state transition
        measurementFcn       % Function handle for the nonlinear measurement
        Q_n_fcn              % Function handle for DMC process noise covariance (free function)
        Q_base               % Base process noise covariance (for other uncertainties)
        R                    % Measurement noise covariance matrix
        x                    % State estimate vector
        P                    % State covariance matrix
        n_orig               % Original state dimension (e.g., position, velocity)
        n_dmc                % DMC state dimension (Gauss-Markov process)
        n                    % Augmented state dimension (n_orig + n_dmc)
        m                    % Measurement dimension
        alpha                % UKF parameter: spread of sigma points (e.g., 1e-3)
        beta_ukf             % UKF parameter: incorporates prior knowledge (e.g., 2)
        kappa                % UKF parameter: secondary scaling (e.g., 0)
        lambda               % UKF scaling parameter
        gamma                % UKF scaling parameter
        Wm                   % Weights for the mean
        Wc                   % Weights for the covariance
        t_prev               % Previous time for time-difference calculation
        dmc_beta             % Beta parameter for the Gauss-Markov process
        dmc_sigma_u2         % Sigma_u^2 parameter for the Gauss-Markov process
    end

    methods
        function obj = UKF_DMC(stateTransitionFcn, measurementFcn, initial_x, initial_P, ...
                               dmc_beta, dmc_sigma_u2, Q_base, R, alpha, beta_ukf, kappa, t_initial)
            %UKF_DMC Constructor for the UKF_DMC class.
            %
            % Inputs:
            %   stateTransitionFcn - Function handle, x_{k+1} = f(x_k, dt).
            %   measurementFcn     - Function handle, z_k     = h(x_k).
            %   initial_x          - Initial augmented state estimate vector.
            %   initial_P          - Initial augmented state covariance matrix.
            %   dmc_beta           - Beta parameter for the Gauss-Markov process.
            %   dmc_sigma_u2       - Sigma_u^2 parameter for the Gauss-Markov process.
            %   Q_base             - Base process covariance matrix (non-DMC part).
            %   R                  - Measurement noise covariance matrix.
            %   alpha, beta_ukf, kappa - UKF tuning parameters.
            %   t_initial          - Initial time for the filter.

            obj.stateTransitionFcn = stateTransitionFcn;
            obj.measurementFcn     = measurementFcn;
            obj.Q_base = Q_base;
            obj.R      = R;
            obj.x      = initial_x;
            obj.P      = initial_P;

            % DMC parameters
            obj.dmc_beta     = dmc_beta;
            obj.dmc_sigma_u2 = dmc_sigma_u2;

            % Set state dimensions (augmented)
            obj.n      = length(initial_x);
            obj.n_orig = obj.n - 3;   % Assuming 3 DMC states are tracked analytically
            obj.n_dmc  = 3;

            % Measurement dimension
            obj.m = size(R, 1);

            % UKF parameters
            obj.alpha   = alpha;
            obj.beta_ukf= beta_ukf;
            obj.kappa   = kappa;

            % Initialize time
            obj.t_prev = t_initial;

            % Pre-calculate weights and scaling
            obj = calculateWeightsAndScaling(obj);

            % Create function handle for the (free) DMC process-noise helper
            % NOTE: this binds beta and sigma_u2 now, and takes only dt at call time.
            obj.Q_n_fcn = @(dt) computeDmcProcessNoise(dt, obj.dmc_beta, obj.dmc_sigma_u2);
        end

        function obj = predict(obj, t_now)
            %PREDICT Performs the prediction step using the unscented transform.
            %
            % Input:
            %   t_now - The current time.

            dt = t_now - obj.t_prev;

            % Generate sigma points from the current state
            sigma_points = generateSigmaPoints(obj);

            % Propagate sigma points through the nonlinear state transition function
            propagated_sigma_points = zeros(obj.n, 2*obj.n + 1);
            for i = 1:(2*obj.n + 1)
                propagated_sigma_points(:, i) = obj.stateTransitionFcn(sigma_points(:, i), dt);
            end

            % Recalculate mean and covariance from propagated sigma points
            obj.x = sum(obj.Wm .* propagated_sigma_points, 2);
            P_pred = zeros(obj.n, obj.n);
            for i = 1:(2*obj.n + 1)
                diff   = propagated_sigma_points(:, i) - obj.x;
                P_pred = P_pred + obj.Wc(i) * (diff * diff.');
            end

            % Build augmented process noise covariance
            Q_augmented = zeros(obj.n, obj.n);

            % Q_base corresponds to original states (e.g., pos, vel, accel)
            Q_augmented(1:3, 1:3) = obj.Q_base;

            % DMC process-noise for the Gauss-Markov acceleration deviation
            Q_dmc = obj.Q_n_fcn(dt);        % 3x3 block (Singer/CA)
            dmc_noise_var = Q_dmc(3, 3);    % relevant term for accel-state channel
            Q_augmented(4, 4) = dmc_noise_var;

            % Update covariance and time
            obj.P     = P_pred + Q_augmented;
            obj.t_prev = t_now;
        end

        function obj = correct(obj, z)
            %CORRECT Performs the correction (update) step of the UKF.

            % Generate sigma points from the predicted state
            sigma_points = generateSigmaPoints(obj);

            % Predict measurements
            Z = zeros(obj.m, 2*obj.n + 1);
            for i = 1:(2*obj.n + 1)
                Z(:, i) = obj.measurementFcn(sigma_points(:, i));
            end

            % Predicted mean measurement and innovation covariances
            z_pred = sum(obj.Wm .* Z, 2);
            Pzz = zeros(obj.m, obj.m);
            Pxz = zeros(obj.n, obj.m);
            for i = 1:(2*obj.n + 1)
                z_diff = Z(:, i) - z_pred;
                x_diff = sigma_points(:, i) - obj.x;
                Pzz = Pzz + obj.Wc(i) * (z_diff * z_diff.');
                Pxz = Pxz + obj.Wc(i) * (x_diff * z_diff.');
            end
            Pzz = Pzz + obj.R;

            % Kalman gain and update
            K     = Pxz / Pzz;
            obj.x = obj.x + K * (z - z_pred);
            obj.P = obj.P - K * Pzz * K.';
        end

        function x_est = getState(obj)
            %GETSTATE Returns the current state estimate.
            x_est = obj.x;
        end

        function P_est = getCovariance(obj)
            %GETCOVARIANCE Returns the current state covariance.
            P_est = obj.P;
        end
    end

    methods (Access = private)
        function obj = calculateWeightsAndScaling(obj)
            %calculateWeightsAndScaling Calculates the unscented transform weights.

            obj.lambda = obj.alpha^2 * (obj.n + obj.kappa) - obj.n;
            obj.gamma  = sqrt(obj.n + obj.lambda);

            obj.Wm = zeros(1, 2*obj.n + 1);
            obj.Wm(1)     = obj.lambda / (obj.n + obj.lambda);
            obj.Wm(2:end) = 1 / (2 * (obj.n + obj.lambda));

            obj.Wc = zeros(1, 2*obj.n + 1);
            obj.Wc(1)     = obj.lambda / (obj.n + obj.lambda) + (1 - obj.alpha^2 + obj.beta_ukf);
            obj.Wc(2:end) = 1 / (2 * (obj.n + obj.lambda));
        end

        function sigma_points = generateSigmaPoints(obj)
            %generateSigmaPoints Generates sigma points from current state mean and covariance.

            % Stabilize P (symmetrize + epsilon on diagonal) before chol
            P_sym = (obj.P + obj.P.') / 2;
            eps_d = 1e-12;
            P_stab = P_sym + eps_d * eye(obj.n);

            A = real(chol(P_stab).');   % upper-triangular, transpose to get columns

            sigma_points = zeros(obj.n, 2*obj.n + 1);
            sigma_points(:, 1) = obj.x;

            for i = 1:obj.n
                sigma_points(:, i + 1)        = obj.x + obj.gamma * A(:, i);
                sigma_points(:, i + obj.n + 1)= obj.x - obj.gamma * A(:, i);
            end
        end
    end
end

% ===== Free helper function (file-local, not a class method) =====
function Q_n = computeDmcProcessNoise(dt, beta, sigma_u2)
%computeDmcProcessNoise Computes the 3×3 DMC (Singer/CA) process noise covariance.
% Handles the constant-acceleration special case when beta ≈ 0.

t = dt;

% Input validation
if sigma_u2 < 0
    error('sigma_u2 must be non-negative.');
end
if t < 0
    error('dt must be non-negative.');
end

if abs(beta) < 1e-9
    % Constant-acceleration limit (beta -> 0)
    Q_n_1_1 = (t^5) / 20;
    Q_n_1_2 = (t^4) / 8;
    Q_n_1_3 = (t^3) / 6;
    Q_n_2_2 = (t^3) / 3;
    Q_n_2_3 = (t^2) / 2;
    Q_n_3_3 = t;
else
    % Singer model (beta > 0)
    e_bt  = exp(-beta * t);
    e_2bt = exp(-2 * beta * t);

    Q_n_1_1 = (1/(3*beta^2))*t^3 + (1/(beta^3))*t^2 + ...
              (1/(beta^4))*(1 - 2*e_bt) + (1/(2*beta^5))*(1 - e_2bt);

    Q_n_1_2 = (1/(2*beta^2))*t^2 + (1/(beta^3))*(1 - e_bt) + ...
              (1/(beta^4))*t*e_bt + (1/(2*beta^5))*(1 - e_2bt);

    Q_n_1_3 = (1/(2*beta))*(1 - e_2bt) + (1/(beta^2))*t*e_bt;

    Q_n_2_2 = (1/(beta^2))*t + (2/(beta^3))*(1 - e_bt) + ...
              (1/(2*beta^3))*(1 - e_2bt);

    Q_n_2_3 = (1/(2*beta^2))*(1 - e_bt)^2;

    Q_n_3_3 = (1/(2*beta))*(1 - e_2bt);
end

% Build symmetric matrix and apply scaling
Q_n = sigma_u2 * [Q_n_1_1, Q_n_1_2, Q_n_1_3;
                  Q_n_1_2, Q_n_2_2, Q_n_2_3;
                  Q_n_1_3, Q_n_2_3, Q_n_3_3];
end