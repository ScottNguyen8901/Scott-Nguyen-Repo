classdef UKF
%UKF  Unscented Kalman Filter (additive-noise) for nonlinear systems.
%
%   This class implements a standard UKF with tuning parameters (alpha,beta,kappa).
%   It assumes ADDITIVE process/measurement noise:
%
%       x_{k+1} = f(x_k, u_k) + w_k,             w_k ~ N(0, Q)
%       z_k     = h(x_k)       + v_k,             v_k ~ N(0, R)
%
%   CONSTRUCTOR
%       obj = UKF(fFcn, hFcn, x0, P0, Q, R, alpha, beta, kappa)
%
%       fFcn   : function handle, state transition, signature f(x[,u])
%                - If you have a control input u, pass it to PREDICT as obj.predict(u)
%       hFcn   : function handle, measurement model, signature h(x)
%       x0     : initial state estimate (n×1)
%       P0     : initial covariance (n×n, PSD)
%       Q, R   : process and measurement noise covariances
%       alpha  : spread of sigma points (1e-3–1e-1 typical)
%       beta   : prior knowledge of distribution (2 for Gaussian optimal)
%       kappa  : secondary scaling (often 0 or 3-n)
%
%   MAIN METHODS
%       obj = predict(obj, u)     % propagate state/covariance through f
%       obj = correct(obj, z)     % incorporate measurement z via h
%       x   = getState(obj)
%       P   = getCovariance(obj)
%
%   NOTES
%       • Uses lower-triangular Cholesky for sigma generation.
%       • Includes mild jitter if P is near-singular (numerical guard).
%       • Cross-covariance uses Wc weights (UKF standard).
%
    properties (Access = private)
        % Models
        fFcn           % state transition: f(x[,u])
        hFcn           % measurement:      h(x)

        % Noise
        Q              % process noise covariance (n×n)
        R              % measurement noise covariance (m×m)

        % Estimates
        x              % current state estimate (n×1)
        P              % current covariance    (n×n)

        % Dimensions
        n              % state dimension
        m              % measurement dimension

        % UKF tuning & derived weights
        alpha
        beta
        kappa
        lambda
        gamma
        Wm             % weights for mean   (1×(2n+1))
        Wc             % weights for cov    (1×(2n+1))
    end

    methods
        function obj = UKF(fFcn, hFcn, x0, P0, Q, R, alpha, beta, kappa)
        %UKF Constructor.
        %   See class header for parameter descriptions.
            arguments
                fFcn (1,1) function_handle
                hFcn (1,1) function_handle
                x0  (:,1) double
                P0  double
                Q   double
                R   double
                alpha (1,1) double = 1e-3
                beta  (1,1) double = 2
                kappa (1,1) double = 0
            end

            % Store models/noise
            obj.fFcn = fFcn;
            obj.hFcn = hFcn;
            obj.x    = x0;
            obj.P    = P0;
            obj.Q    = Q;
            obj.R    = R;

            % Dimensions
            obj.n = numel(x0);
            obj.m = size(R,1);

            % Tuning
            obj.alpha = alpha;
            obj.beta  = beta;
            obj.kappa = kappa;

            % Precompute weights/scales
            obj = obj.computeWeights();
        end

        function obj = predict(obj, u)
        %PREDICT Propagate sigma points through f and update (x,P) with Q.
        %   obj = predict(obj)       % no control
        %   obj = predict(obj, u)    % with control input u (passed to f)
            if nargin < 2, u = []; end

            % Generate sigma points from current (x,P)
            Xi = obj.sigmaPoints(obj.x, obj.P);   % n×(2n+1)

            % Propagate through dynamics
            Yi = zeros(obj.n, size(Xi,2));
            if isempty(u)
                for i = 1:size(Xi,2)
                    Yi(:,i) = obj.fFcn(Xi(:,i));
                end
            else
                for i = 1:size(Xi,2)
                    Yi(:,i) = obj.fFcn(Xi(:,i), u);
                end
            end

            % Recompute mean/cov
            x_pred = Yi * obj.Wm.';                           % n×1
            P_pred = zeros(obj.n);
            for i = 1:size(Yi,2)
                d = Yi(:,i) - x_pred;
                P_pred = P_pred + obj.Wc(i) * (d * d.');
            end
            P_pred = P_pred + obj.Q;                          % add process noise

            % Commit
            obj.x = x_pred;
            obj.P = obj.symmetrize(P_pred);
        end

        function obj = correct(obj, z)
        %CORRECT Update (x,P) using measurement z through h and R.
            % Sigma points from predicted (x,P)
            Xi = obj.sigmaPoints(obj.x, obj.P);               % n×(2n+1)

            % Predict measurements
            Zi = zeros(obj.m, size(Xi,2));
            for i = 1:size(Xi,2)
                Zi(:,i) = obj.hFcn(Xi(:,i));
            end
            z_pred = Zi * obj.Wm.';                           % m×1

            % Innovation covariance S and cross-covariance Pxz
            S   = zeros(obj.m);
            Pxz = zeros(obj.n, obj.m);
            for i = 1:size(Xi,2)
                dx = Xi(:,i) - obj.x;
                dz = Zi(:,i) - z_pred;
                S   = S   + obj.Wc(i) * (dz * dz.');
                Pxz = Pxz + obj.Wc(i) * (dx * dz.');
            end
            S = S + obj.R;

            % Kalman gain (solve instead of invert for stability)
            % K = Pxz / S
            K = Pxz / (S + 1e-12*eye(obj.m));

            % Update
            obj.x = obj.x + K * (z - z_pred);
            obj.P = obj.symmetrize(obj.P - K * S * K.');
        end

        function x = getState(obj)
        %GETSTATE Current state estimate (n×1).
            x = obj.x;
        end

        function P = getCovariance(obj)
        %GETCOVARIANCE Current state covariance (n×n).
            P = obj.P;
        end
    end

    methods (Access = private)
        function obj = computeWeights(obj)
        %COMPUTEWEIGHTS Pre-compute lambda, gamma, and UKF weights.
            obj.lambda = obj.alpha^2 * (obj.n + obj.kappa) - obj.n;
            obj.gamma  = sqrt(obj.n + obj.lambda);

            obj.Wm = zeros(1, 2*obj.n + 1);
            obj.Wc = zeros(1, 2*obj.n + 1);

            obj.Wm(1) = obj.lambda / (obj.n + obj.lambda);
            obj.Wc(1) = obj.Wm(1) + (1 - obj.alpha^2 + obj.beta);

            obj.Wm(2:end) = 1/(2*(obj.n + obj.lambda));
            obj.Wc(2:end) = obj.Wm(2:end);
        end

        function Xi = sigmaPoints(obj, x, P)
        %SIGMAPOINTS Generate sigma points around mean x with covariance P.
        %   Returns n×(2n+1) matrix: [x, x±gamma*A_i], A=chol(P,'lower')
            % Ensure P is numerically SPD
            jitter = 1e-12;
            Aok = false; tries = 0;
            while ~Aok && tries < 5
                try
                    A = chol(P + jitter*eye(obj.n), 'lower');
                    Aok = true;
                catch
                    jitter = max(1e-12, 10*jitter);
                    tries  = tries + 1;
                end
            end
            if ~Aok
                % final fallback: make symmetric PSD
                [V,D] = eig((P+P')/2);
                D = max(D, 1e-12*eye(obj.n));
                A = chol(V*D*V','lower');
            end

            Xi = zeros(obj.n, 2*obj.n + 1);
            Xi(:,1) = x;
            for i = 1:obj.n
                Xi(:,1+i)      = x + obj.gamma * A(:,i);
                Xi(:,1+obj.n+i)= x - obj.gamma * A(:,i);
            end
        end

        function M = symmetrize(~, M)
        %SYMMETRIZE Enforce symmetry (numerical hygiene).
            M = (M + M.')/2;
        end
    end
end