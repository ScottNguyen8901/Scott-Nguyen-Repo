function [q, R_bi_QUEST, J] = quest(v_b, v_i, w_k)
    % DESCRIPTION
    %   Computes the optimal quaternion, rotation matrix, and Wahba cost
    %   function using the QUEST algorithm.
    %
    % INPUTS       size         Type        Description
    %   v_b        (3, n)       Double      Vectors in the body frame
    %   v_i        (3, n)       Double      Vectors in the inertial frame
    %   w_k        (n, 1)       Double      Weights for each vector pair
    %
    % OUTPUTS      size         Type        Description
    %   q          (4, 1)       Double      Optimal quaternion in the form [qx, qy, qz, qw]'
    %   R_bi_QUEST (3, 3)       Double      Optimal rotation matrix from body to inertial frame
    %   J          (1, 1)       Double      Wahba's cost function value
    %
    % FUNCTION
    
    % Number of vector pairs
    n = length(w_k);
    
    % Initialize B matrix
    B = zeros(3);
    
    % Compute the weighted outer product sum
    for i = 1:n
        B = B + w_k(i) * (v_b(:,i) * v_i(:,i)');  % Outer product of columns
    end
    
    % Symmetric part of B
    S = B + B';
    Z = [B(2,3) - B(3,2);
         B(3,1) - B(1,3);
         B(1,2) - B(2,1)];
    sigma = trace(B);
    
    % Compute lambda_opt as the sum of weights
    lambda_opt = sum(w_k);
    
    % Solve for vector p
    p = ((lambda_opt + sigma) * eye(3) - S) \ Z;
    
    % Compute the quaternion q
    q = (1 / sqrt(1 + p' * p)) * [p; 1];
    
    % Convert quaternion to rotation matrix
    R_bi_QUEST = quat_to_rot(q);
    
    % Compute Wahba's cost function
    J = wahba_cost(v_b, v_i, R_bi_QUEST, w_k);
end
