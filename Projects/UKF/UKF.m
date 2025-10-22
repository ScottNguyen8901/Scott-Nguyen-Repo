classdef UKF
    %UKF A custom class for implementing an unscented Kalman filter.
    % This implementation assumes additive process and measurement noise
    % and requires user-defined state transition and measurement functions.

    properties (Access = private)
        stateTransitionFcn   % Function handle for the nonlinear state transition
        measurementFcn       % Function handle for the nonlinear measurement
        Q                    % Process noise covariance matrix
        R                    % Measurement noise covariance matrix
        x                    % State estimate vector
        P                    % State covariance matrix
        n                    % State dimension
        m                    % Measurement dimension
        alpha                % UKF parameter: spread of sigma points (e.g., 1e-3)
        beta                 % UKF parameter: incorporates prior knowledge (e.g., 2)
        kappa                % UKF parameter: secondary scaling (e.g., 0)
        lambda               % UKF scaling parameter
        gamma                % UKF scaling parameter
        Wm                   % Weights for the mean
        Wc                   % Weights for the covariance
    end

    methods
        function obj = UKF(stateTransitionFcn, measurementFcn, initial_x, initial_P, Q, R, alpha, beta, kappa)
            %UKF Constructor for the UKF class.
            % Initializes the filter with nonlinear functions and initial conditions.
            %
            % Inputs:
            %   stateTransitionFcn - Function handle,   x_{k+1} = f(x_k).
            %   measurementFcn     - Function handle,   z_k     = h(x_k).
            %   initial_x          - Initial state estimate vector.
            %   initial_P          - Initial state covariance matrix.
            %   Q                  - Process noise covariance matrix.
            %   R                  - Measurement noise covariance matrix.
            %   alpha, beta, kappa - UKF tuning parameters.

            obj.stateTransitionFcn = stateTransitionFcn;
            obj.measurementFcn     = measurementFcn;
            obj.Q = Q;
            obj.R = R;
            obj.x = initial_x;
            obj.P = initial_P;

            % State dimension
            obj.n = length(initial_x);
            % Measurement dimension
            obj.m = size(R, 1);

            % Set UKF parameters
            obj.alpha = alpha;
            obj.beta  = beta;
            obj.kappa = kappa;

            % Pre-calculate UKF weights and scaling factors
            obj = calculateWeightsAndScaling(obj);
        end

        function obj = predict(obj)
            %PREDICT Performs the prediction step using the unscented transform.

            % Generate sigma points
            sigma_points = generateSigmaPoints(obj);

            % Propagate sigma points through the nonlinear state transition function
            propagated_sigma_points = zeros(obj.n, 2*obj.n + 1);
            for i = 1:(2*obj.n + 1)
                propagated_sigma_points(:, i) = obj.stateTransitionFcn(sigma_points(:, i));
            end

            % Recalculate mean and covariance from propagated sigma points
            obj.x = sum(obj.Wm .* propagated_sigma_points, 2);
            P_pred = zeros(obj.n, obj.n);
            for i = 1:(2*obj.n + 1)
                diff   = propagated_sigma_points(:, i) - obj.x;
                P_pred = P_pred + obj.Wc(i) * (diff * diff.');
            end
            obj.P = P_pred + obj.Q;  % Add process noise covariance
        end

        function obj = correct(obj, z)
            %CORRECT Performs the correction (update) step of the UKF.
            %
            % Inputs:
            %   z  - The new measurement vector.

            % Generate sigma points from the predicted state
            sigma_points = generateSigmaPoints(obj);

            % Predict measurements by propagating sigma points through measurement function
            predicted_measurements = zeros(obj.m, 2*obj.n + 1);
            for i = 1:(2*obj.n + 1)
                predicted_measurements(:, i) = obj.measurementFcn(sigma_points(:, i));
            end

            % Recalculate predicted mean measurement and innovation covariance
            z_pred = sum(obj.Wm .* predicted_measurements, 2);
            Pzz = zeros(obj.m, obj.m);
            Pxz = zeros(obj.n, obj.m);
            for i = 1:(2*obj.n + 1)
                z_diff = predicted_measurements(:, i) - z_pred;
                x_diff = sigma_points(:, i) - obj.x;
                Pzz = Pzz + obj.Wc(i) * (z_diff * z_diff.');
                Pxz = Pxz + obj.Wc(i) * (x_diff * z_diff.');
            end
            Pzz = Pzz + obj.R;  % Add measurement noise covariance

            % Calculate Kalman gain and update state and covariance
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

            % Calculate scaling parameters
            obj.lambda = obj.alpha^2 * (obj.n + obj.kappa) - obj.n;
            obj.gamma  = sqrt(obj.n + obj.lambda);

            % Calculate weights for mean and covariance
            obj.Wm = zeros(1, 2*obj.n + 1);
            obj.Wm(1)    = obj.lambda / (obj.n + obj.lambda);
            obj.Wm(2:end)= 1 / (2 * (obj.n + obj.lambda));

            obj.Wc = zeros(1, 2*obj.n + 1);
            obj.Wc(1)    = obj.lambda / (obj.n + obj.lambda) + (1 - obj.alpha^2 + obj.beta);
            obj.Wc(2:end)= 1 / (2 * (obj.n + obj.lambda));
        end

        function sigma_points = generateSigmaPoints(obj)
            %generateSigmaPoints Generates sigma points from current state mean and covariance.

            A = real(chol(obj.P).');  % Cholesky decomposition of covariance (upper)

            % Allocate
            sigma_points = zeros(obj.n, 2*obj.n + 1);

            % Sigma point 0 (mean)
            sigma_points(:, 1) = obj.x;

            % Sigma points 1 to n (positive scaled deviations)
            for i = 1:obj.n
                sigma_points(:, i + 1) = obj.x + obj.gamma * A(:, i);
            end

            % Sigma points n+1 to 2n (negative scaled deviations)
            for i = 1:obj.n
                sigma_points(:, i + obj.n + 1) = obj.x - obj.gamma * A(:, i);
            end
        end
    end
end