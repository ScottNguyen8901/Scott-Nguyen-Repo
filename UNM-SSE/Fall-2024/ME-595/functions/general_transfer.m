function [dv_1, dv_2, dv_t] = general_transfer(r_1, r_2, p, e, mu)
%
% DESCRIPTION
%   Calculate the delta-Vs and total delta-V for a general orbital transfer 
%   between two circular orbits, considering an elliptical transfer orbit. 
%   This function allows for a specified semi-latus rectum and eccentricity 
%   of the transfer orbit, which is useful for non-Hohmann transfers.
%
% INPUTS        size    Type
%   r_1         (1,1)   Double  Radius of initial circular orbit [DU]
%   r_2         (1,1)   Double  Radius of final circular orbit   [DU]
%   p           (1,1)   Double  Semi-latus rectum of transfer orbit [DU]
%   e           (1,1)   Double  Eccentricity of transfer orbit
%   mu          (1,1)   Double  Gravitational parameter [DU^3/TU^2]
%
% OUTPUTS       size    Type
%   dv_1        (1,1)   Double  First delta-V  [DU/TU]
%   dv_2        (1,1)   Double  Second delta-V [DU/TU]
%   dv_t        (1,1)   Double  Total delta-V  [DU/TU]
%
% FUNCTION
%

    E_t = -mu * (1 - e^2) / (2 * p);  % Energy of the transfer orbit
    h_t = sqrt(mu * p);               % Angular momentum of the transfer orbit
    
    % First delta-V at the initial orbit
    nu_1 = acosd((p - r_1) / (r_1 * e));  % True anomaly at the initial orbit
    v_c_1 = sqrt(mu / r_1);               % Circular orbit velocity at initial orbit
    v_1 = sqrt(2 * (E_t + mu / r_1));     % Transfer orbit velocity at r_1
    v_1_cos_phi_1 = h_t / r_1;
    dv_1 = sqrt(v_1^2 + v_c_1^2 - 2 * v_c_1 * v_1_cos_phi_1);
    
    % Second delta-V at the final orbit
    nu_2 = acosd((p - r_2) / (r_2 * e));  % True anomaly at the final orbit
    v_c_2 = sqrt(mu / r_2);               % Circular orbit velocity at final orbit
    v_2 = sqrt(2 * (E_t + mu / r_2));     % Transfer orbit velocity at r_2
    v_2_cos_phi_2 = h_t / r_2;
    dv_2 = sqrt(v_2^2 + v_c_2^2 - 2 * v_c_2 * v_2_cos_phi_2);
    
    % Total delta-V for the transfer
    dv_t = dv_1 + dv_2;

end