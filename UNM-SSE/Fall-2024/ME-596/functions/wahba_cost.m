function J = wahba_cost(v_b, v_i, R_bi_q, w_k)
    % DESCRIPTION
    %   Computes a weighted sum of squared differences between vectors in 
    %   the body frame and their rotated counterparts in the inertial frame.
    %
    % INPUTS       size         Type        Description
    %   v_b        (3, n)       Double      Vectors in the body frame
    %   v_i        (3, n)       Double      Vectors in the inertial frame
    %   R_bi_q     (3, 3)       Double      Rotation matrix from body to inertial frame
    %   w_k        (n, 1)       Double      Weights for each vector pair
    %   n          (1, 1)       Integer     Number of vector pairs
    %
    % OUTPUTS      size         Type        Description
    %   J          (1, 1)       Double      Weighted sum of squared differences
    %
    % FUNCTION
    
    n = length(w_k);

    % Initialize J
    J = 0;
    
    % Loop over each vector pair and accumulate the weighted squared differences
    for i = 1:n
        diff = v_b(:,i) - R_bi_q * v_i(:,i);  % Difference between the vectors
        J = J + w_k(i) * (diff' * diff);      % Weighted sum of squared differences
    end
end
