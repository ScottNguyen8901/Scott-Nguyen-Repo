function R = rot_mat(W, w, i)
    %
    % DESCRIPTION
    %   Computes the rotation matrix for transforming from the perifocal frame 
    %   to the ECI (Earth-Centered Inertial) or heliocentric frame based on the orbital 
    %   elements: Right Ascension of the Ascending Node (W), Argument of Periapsis (w), 
    %   and Inclination (i).
    %   This function assumes that the orbital elements are given relative to the 
    %   Earth-centered equatorial coordinate system or the ecliptic plane for heliocentric orbits.
    %
    % INPUTS         size         Type     Description                               Units
    %   W            (1,1)        Double   Right Ascension of Ascending Node         [rad]
    %   w            (1,1)        Double   Argument of Periapsis                    [rad]
    %   i            (1,1)        Double   Inclination                              [rad]
    %
    % OUTPUTS        size         Type     Description                               Units
    %   R            (3,3)        Double   3x3 rotation matrix from the perifocal frame
    %                                         to the ECI (or heliocentric) frame
    %
    % NOTES
    %   The function computes the rotation matrix from the perifocal frame to the ECI frame
    %   or heliocentric frame based on the orbital elements. It uses the appropriate rotation
    %   matrices for the transformations: first for the RAAN (W), then for the inclination (i),
    %   and finally for the argument of periapsis (w).
    %   The returned matrix is used for converting position and velocity vectors from the 
    %   perifocal coordinate system to the ECI or heliocentric frame.
    %

    % Rotation around the Z-axis by W (RAAN)
    R1 = [cos(W), sin(W), 0;
          -sin(W), cos(W), 0;
          0, 0, 1];

    % Rotation around the X-axis by i (Inclination)
    R2 = [1, 0, 0;
          0, cos(i), sin(i);
          0, -sin(i), cos(i)];

    % Rotation around the Z-axis by w (Argument of Periapsis)
    R3 = [cos(w), sin(w), 0;
         -sin(w), cos(w), 0;
          0, 0, 1];

    % Combined rotation matrix
    R = R1 * R2 * R3;
end
