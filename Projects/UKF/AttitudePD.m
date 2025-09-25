classdef AttitudePD
%AttitudePD  Simple body-frame PD attitude controller (quaternion-based).
%
%   Usage (typical):
%       J  = diag([10 8 5]);             % inertia (kg·m^2) if you want auto-tune
%       ctl = AttitudePD('AutoTune',true,'J',J, ...
%                        'Zeta',0.9,'OmegaN',0.8);  % optional auto-tuning
%       % or manual gains:
%       % ctl = AttitudePD('Kp',diag([8 8 6]), 'Kd',diag([6 6 5]));
%
%       ctl = ctl.setDesired([1;0;0;0]);  % desired quaternion (scalar-first)
%       tau = ctl.torque(q_current, w_current);
%
%   Torque law:
%       tau = -Kp * q_ev - Kd * w
%     where q_e = q_d^{-1} ⊗ q  (shortest-rotation, i.e., if q_e0<0, flip sign)
%
%   Notes:
%     • Quaternion convention: scalar-first [q0;q1;q2;q3], unit-norm.
%     • All vectors are body-frame; ω in rad/s; τ in N·m.
%     • You may enable saturation via TauMax (per-axis).
%
%   Parameters (Name-Value pairs in constructor):
%       'Kp'      : 3x3, 3x1, or scalar proportional gain (default: auto if AutoTune)
%       'Kd'      : 3x3, 3x1, or scalar derivative gain (default: auto if AutoTune)
%       'AutoTune': true/false (default: false). If true, requires J and sets
%                   Kp_i = J_i * OmegaN^2,  Kd_i = 2*Zeta*J_i*OmegaN.
%       'J'       : 3x3 inertia (needed when AutoTune=true)
%       'OmegaN'  : desired natural frequency (rad/s), default 0.8
%       'Zeta'    : damping ratio, default 0.9
%       'TauMax'  : 3x1 nonnegative saturation limits (N·m). Empty = no sat.
%
%   Methods:
%       setDesired(qd)   : set desired quaternion (scalar-first)
%       torque(q, w)     : compute τ given current attitude & rate
%       setGains(Kp,Kd)  : update gains
%
%   (c) you — keep it simple :)

    properties
        Kp   double {mustBeNumeric} = diag([8 8 8])     % 3x3 (or scalar/3x1 accepted)
        Kd   double {mustBeNumeric} = diag([6 6 6])     % 3x3 (or scalar/3x1 accepted)
        qd   double {mustBeVector, mustBeFinite} = [1;0;0;0]   % desired quaternion (unit)
        TauMax double = []                               % 3x1 per-axis saturation (optional)
    end

    properties (Access=private)
        useSaturation logical = false
    end

    methods
        function obj = AttitudePD(varargin)
            % Parse name-value args
            p = inputParser;
            addParameter(p,'Kp',[]);
            addParameter(p,'Kd',[]);
            addParameter(p,'AutoTune',false,@islogical);
            addParameter(p,'J',[]);
            addParameter(p,'OmegaN',0.8,@(x)isscalar(x)&&x>0);
            addParameter(p,'Zeta',0.9,@(x)isscalar(x)&&x>0);
            addParameter(p,'TauMax',[]);
            parse(p,varargin{:});
            S = p.Results;

            % Auto-tune if requested
            if S.AutoTune
                assert(~isempty(S.J), 'AutoTune requires inertia matrix J.');
                J = diag(S.J);  % assume principal axes (diagonal J)
                wn = S.OmegaN;  z = S.Zeta;
                kp = J * wn^2;
                kd = 2*z*J*wn;
                obj.Kp = diag(kp);
                obj.Kd = diag(kd);
            else
                if ~isempty(S.Kp), obj.Kp = AttitudePD.to3x3(S.Kp); end
                if ~isempty(S.Kd), obj.Kd = AttitudePD.to3x3(S.Kd); end
            end

            % Optional saturation
            if ~isempty(S.TauMax)
                obj = obj.setSaturation(S.TauMax);
            end
        end

        function obj = setDesired(obj, qd)
            obj.qd = qd(:) / norm(qd);
        end

        function obj = setGains(obj, Kp, Kd)
            obj.Kp = AttitudePD.to3x3(Kp);
            obj.Kd = AttitudePD.to3x3(Kd);
        end

        function obj = setSaturation(obj, tauMax)
            if isempty(tauMax)
                obj.useSaturation = false; obj.TauMax = [];
            else
                tauMax = abs(tauMax(:));
                assert(numel(tauMax)==3,'TauMax must be 3x1 or empty.');
                obj.useSaturation = true; obj.TauMax = tauMax;
            end
        end

        function tau = torque(obj, q, w)
            %TORQUE  Compute body torque τ = -Kp*q_ev - Kd*w toward desired qd.
            % Inputs:
            %   q : current quaternion (scalar-first, unit)
            %   w : current body rates (rad/s)
            % Output:
            %   tau : 3x1 body torque (N·m)
            q = q(:) / norm(q);
            w = w(:);

            % Error quaternion q_e = qd^{-1} ⊗ q  (qd_inv = [qd0;-qdv])
            qd0 = obj.qd(1); qdv = obj.qd(2:4);
            qinv = [qd0; -qdv];

            % Quaternion multiply q_e = qinv ⊗ q (scalar-first)
            s1=qinv(1); v1=qinv(2:4);
            s2=q(1);    v2=q(2:4);
            q_e0 = s1*s2 - dot(v1,v2);
            q_ev = s1*v2 + s2*v1 + cross(v1,v2);

            % Shortest-rotation convention (avoid unwinding)
            if q_e0 < 0
                q_e0 = -q_e0; 
                q_ev = -q_ev;
            end

            % PD torque
            tau = - (obj.Kp * q_ev + obj.Kd * w);

            % Optional saturation
            if obj.useSaturation
                tau = max(min(tau, obj.TauMax), -obj.TauMax);
            end
        end
    end

    methods (Static, Access=private)
        function M = to3x3(G)
            % Accept scalar, 3x1, or 3x3 gains
            if isscalar(G)
                M = eye(3)*G;
            elseif isvector(G) && numel(G)==3
                M = diag(G(:));
            else
                assert(all(size(G)==[3 3]),'Gain must be scalar, 3x1, or 3x3.');
                M = G;
            end
        end
    end
end