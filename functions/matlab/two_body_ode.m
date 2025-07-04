function dydt = two_body_ode(t, y)
%
% DESCRIPTION
%   Computes the derivative for two-body motion
%
% INPUTS   Size     Type      Description                       Units
%   t      (1,1)    (Double)   Time                             [TU]
%   y      (6,1)    (Double)   State vector: [pos; vel]         [DU, DU/TU]
%
% OUTPUTS  Size   Type      Description                    Units
%   dydt   (6,1)  (Double)  Derivative of the state vector:[DU/TU, DU/TU^2]
%
% FUNCTION

    constants;
    
    r = y(1:3);
    v = y(4:6);
    
    r_norm = norm(r);
    a = -mu_E * r / r_norm^3;
    dydt = [v; a];
end