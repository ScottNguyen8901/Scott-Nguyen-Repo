function [q, R_bi_q, J, eig_max] = q_method(v_b, v_i, w_k)
    % DESCRIPTION
    %   Compute the quaternion q that transforms vectors from the body frame (b) 
    %   to the inertial frame (i) using the q-method. The function first computes 
    %   the weighted outer product of the vectors and then calculates the quaternion 
    %   that represents the rotation between the frames.
    %
    %   The function also computes the corresponding rotation matrix (R_bi_q) and
    %   a cost function (J) that represents the alignment between the transformed vectors.
    %
    % INPUTS        size    Type
    %   v_b        (3, N)  Double  Matrix of vectors in the body frame
    %   v_i        (3, N)  Double  Matrix of vectors in the inertial frame
    %   w_k        (1, N)  Double  Vector of weights (default: ones if not provided)
    %
    % OUTPUTS       size    Type
    %   q          (4, 1)  Double  Quaternion representing the rotation from body to inertial frame
    %   R_bi_q     (3, 3)  Double  Rotation matrix corresponding to the quaternion q
    %   J          (1, 1)  Double  Cost function representing the alignment of transformed vectors
    %   eig_max    (1, 1)  Double  Maximum eigen value from K matrix
    %
    % FUNCTION

    if nargin < 3
        w_k = ones(1, size(v_b, 2));  % Default weights to 1 if not provided
    end
    
    [~, n] = size(v_b);
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

    % Construct matrix K and compute eigenvalues/eigenvectors
    K = [S - sigma * eye(3), Z;
         Z', sigma];
    [D, V] = eig(K);

    % Eigenvector corresponding to the largest eigenvalue
    [eig_max, idx] = max(diag(V));
    q = D(:, idx);
    
    R_bi_q = quat_to_rot(q);
    
    J = wahba_cost(v_b, v_i, R_bi_q, w_k);

end