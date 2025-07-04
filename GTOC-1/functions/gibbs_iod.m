function [v_1_vec, v_2_vec, v_3_vec] = gibbs_iod(r_1_vec, r_2_vec, r_3_vec, mu)
    %
    % DESCRIPTION
    %   Compute the velocity vectors at three different points using Gibbs' method for orbit determination.
    %
    % INPUTS    size    Type
    %   r_1_vec  (3,1)   Double  Position vector at time 1 [DU]
    %   r_2_vec  (3,1)   Double  Position vector at time 2 [DU]
    %   r_3_vec  (3,1)   Double  Position vector at time 3 [DU]
    %   mu       (1,1)   Double  Gravitational parameter [DU^3/TU^2]
    %
    % OUTPUTS   size    Type
    %   v_1_vec  (3,1)   Double  Velocity vector at time 1 [DU/TU]
    %   v_2_vec  (3,1)   Double  Velocity vector at time 2 [DU/TU]
    %   v_3_vec  (3,1)   Double  Velocity vector at time 3 [DU/TU]
    %
    % FUNCTION
    r_1 = norm(r_1_vec);
    r_2 = norm(r_2_vec);
    r_3 = norm(r_3_vec);
    
    D_vec = cross(r_1_vec, r_2_vec) + ...
            cross(r_2_vec, r_3_vec) + ...
            cross(r_3_vec, r_1_vec);
    
    N_vec = r_3*cross(r_1_vec, r_2_vec) + ...
            r_1*cross(r_2_vec, r_3_vec) + ...
            r_2*cross(r_3_vec, r_1_vec);
    
    S_vec = (r_2 - r_3)*r_1_vec +...
            (r_3 - r_1)*r_2_vec +...
            (r_1 - r_2)*r_3_vec;
    
    D = norm(D_vec);
    N = norm(N_vec);
    S = norm(S_vec);
    
    p = N/D;
    e = S/D;
    
    Q_hat = S_vec/S;
    W_hat = N_vec/N;
    P_hat = cross(Q_hat, W_hat);
    
    L = sqrt(mu / (D * N));
    
    B_1_vec = cross(D_vec, r_1_vec);
    B_2_vec = cross(D_vec, r_2_vec);
    B_3_vec = cross(D_vec, r_3_vec);
    
    v_1_vec = (L / r_1) * B_1_vec + L * S_vec;
    v_2_vec = (L / r_2) * B_2_vec + L * S_vec;
    v_3_vec = (L / r_3) * B_3_vec + L * S_vec;
end