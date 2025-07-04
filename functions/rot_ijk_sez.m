function R_IJK_SEZ = rot_ijk_sez(L, theta)
    %
    % DESCRIPTION
    %   Calculate rotation matrix from IJK to SEZ frame
    %
    % INPUTS       size     Type
    %   L          (1,1)    Double  Longitude [rad]
    %   theta      (1,1)    Double  Latitude  [rad]
    %
    % OUTPUTS
    %   R_IJK_SEZ  (1,1)   Rotation matrix from IJK to SEZ []
    %
    % NOTES
    %   The node line vector is assumed to be along the x-axis for this calculation
    %
    % FUNCTION 

    R_IJK_SEZ = [sin(L)*cos(theta), sin(L)*sin(theta), -cos(L);
                 -sin(theta), cos(theta), 0;
                 cos(L)*cos(theta), cos(L)*sin(theta), sin(L)];
end
