function [R_bi_tri, J] = triad(v_b, v_i, choice)
    % DESCRIPTION
    %   Compute the rotation matrix R_bi that transforms vectors from the
    %   body frame (b) to the inertial frame (i). The function first 
    %   computes orthonormal basis vectors for both frames and then calculates 
    %   the rotation matrix that relates these two frames.
    %
    % INPUTS        size    Type
    %   vb         (3,2)   Double  Matrix containing two vectors in the body frame
    %   vi         (3,2)   Double  Matrix containing two vectors in the inertial frame
    %   choice     (1)     Integer  1 or 2 to control vector order
    %
    % OUTPUTS       size    Type
    %   R_bi       (3,3)   Double  Rotation matrix from body frame to inertial frame
    %
    % FUNCTION
    
    % Check the choice and adjust the vectors' order accordingly
    
    w_k = [1; 1];

    if choice == 1
        % v_1b is the first column, v_2b is the second column
        v_1b = v_b(:,1);
        v_2b = v_b(:,2);
        
        v_1i = v_i(:,1);
        v_2i = v_i(:,2);
    elseif choice == 2
        % Swap the order of v_1b and v_2b, and v_1i and v_2i
        v_1b = v_b(:,2);
        v_2b = v_b(:,1);
        
        v_1i = v_i(:,2);
        v_2i = v_i(:,1);
    else
        error('Choice must be 1 or 2.');
    end

    % Body frame orthonormal vectors
    t_1b = v_1b;
    t_2b = cross(v_1b, v_2b) / norm(cross(v_1b, v_2b));
    t_3b = cross(t_1b, t_2b);

    % Inertial frame orthonormal vectors
    t_1i = v_1i;
    t_2i = cross(v_1i, v_2i) / norm(cross(v_1i, v_2i));
    t_3i = cross(t_1i, t_2i);

    % Rotation matrix from body frame to inertial frame
    R_bt = [t_1b, t_2b, t_3b];
    R_ti = [t_1i, t_2i, t_3i]';

    % Final rotation matrix
    R_bi_tri = R_bt * R_ti;
    
    J = wahba_cost(v_b, v_i, R_bi_tri, w_k);
end
