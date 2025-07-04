function [axes_new, mois_new] = reorder_axes_and_mois(I)
    % DESCRIPTION
    %   Reorders the eigenvectors (axes) and eigenvalues (moments of inertia) 
    %   based on alignment with [1; 0; 0], [0; 1; 0], and [0; 0; 1].
    %
    % INPUTS        size    Type    Description
    %   I           (3,3)   Double  Inertia matrix of the spacecraft [kg*m^2]
    %
    % OUTPUTS          size    Type    Description
    %   AXES_reordered (3,3)   Double  Reordered eigenvectors
    %   MOIS_reordered (3,1)   Double  Reordered eigenvalues
    %
    % FUNCTION
    %   Reorders the principal axes of inertia and the corresponding moments
    %   of inertia such that the first axis aligns closest to [1; 0; 0], the
    %   second to [0; 1; 0], and the third to [0; 0; 1].
    
    % Computing eigenvector and eigenvalues
    [axes, mois] = eig(I);

    % Desired axes for comparison
    desired_axes = [1, 0, 0;   % Closest to [1; 0; 0]
                    0, 1, 0;   % Closest to [0; 1; 0]
                    0, 0, 1];  % Closest to [0; 0; 1]

    % Compute alignments between current axes and desired axes
    alignments = abs(desired_axes * axes); % Take absolute values for alignment

    % Find the best match for each desired axis
    [~, new_order] = max(alignments, [], 2);

    % Reorder the axes and eigenvalues
    axes_new = axes(:, new_order);   % Reorder the eigenvectors
    MOIS_diag = diag(mois);          % Convert MOIS matrix to vector
    mois_new = MOIS_diag(new_order); % Reorder eigenvalues
    mois_new = diag(mois_new);
    
    % Ensure all values in AXES_reordered are positive
    for i = 1:3
        if any(axes_new(i, i) < 0)
            axes_new(:, i) = -axes_new(:, i); % Flip column if negative
        end
    end

end
