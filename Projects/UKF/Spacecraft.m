classdef Spacecraft
%Spacecraft  Rigid-body attitude model (no sim/plot here).
%   Dynamics in body frame:
%       J * wdot = tau(t,x) - w × (J w)
%       qdot     = 0.5 * Omega(w) * q
%   Quaternion convention: scalar-first q = [q0;q1;q2;q3].

    properties
        J               % 3x3 inertia matrix (kg·m^2)
        Jinv            % precomputed inverse of J
        torqueFcn       % handle @(t,x) -> 3x1 tau in body frame (N·m)
    end

    methods
        function obj = Spacecraft(J, torque)
        % Constructor
        %   J      : 3x3 inertia matrix
        %   torque : either 3x1 constant vector OR @(t,x) -> 3x1
            arguments
                J (3,3) double
                torque {mustBeVectorOrFcn}
            end
            obj.J    = J;
            obj.Jinv = inv(J);

            if isa(torque,'function_handle')
                obj.torqueFcn = torque;
            else
                tau = torque(:);
                obj.torqueFcn = @(t,x) tau; %#ok<ASGLU>
            end
        end

        function dx = attitudeDynamics(obj, t, x)
        %ATTITUDEDYNAMICS  State derivative for [w; q].
        %   x = [w(1:3); q(1:4)] with scalar-first quaternion.
            w = x(1:3);
            q = x(4:7);

            tau  = obj.torqueFcn(t, x);                 % N·m
            wdot = obj.Jinv * (tau - cross(w, obj.J*w));

            % qdot = 0.5 * Omega(w) * q  (scalar-first)
            Om = [ 0,    -w.';
                   w,  -skew(w) ];
            qdot = 0.5 * Om * q;

            dx = [wdot; qdot];
        end

        function qn = normalizeQuaternion(~, q)
        %NORMALIZEQUATERNION  Utility to unit-normalize a quaternion.
            n  = norm(q);
            if n == 0, qn = [1;0;0;0]; else, qn = q / n; end
        end
    end
end

% ---------- helpers ----------
function S = skew(v)
    S = [   0, -v(3),  v(2);
          v(3),    0, -v(1);
         -v(2), v(1),    0 ];
end

function mustBeVectorOrFcn(a)
    if ~(isa(a,'function_handle') || (isvector(a) && numel(a)==3))
        error('torque must be a 3x1 vector or a function handle @(t,x).');
    end
end