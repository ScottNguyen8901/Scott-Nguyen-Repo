function R_ECI_SEZ = rot_eci_sez(L, theta)
%
% DESCRIPTION
%   Calculate rotation matrix from IJK to SEZ frame
%
% INPUTS       size     Type     Description                      Units
%   L          (1,1)    Double   Longitude                        [rad]
%   theta      (1,1)    Double   Latitude                         [rad]
%
% OUTPUTS      size     Type     Description                        Units
%   R_IJK_SEZ  (3,3)    Double   Rotation matrix from IJK to SEZ    []
%
% FUNCTION

    R_ECI_SEZ = [sin(L)*cos(theta), sin(L)*sin(theta), -cos(L);
                 -sin(theta), cos(theta), 0;
                 cos(L)*cos(theta), cos(L)*sin(theta), sin(L)];
end