function [FOV_W, FOV_H] = calc_FOV(f, cam_W, cam_H, pix_W, pix_H)
%
% DESCRIPTION
%   Calculates the horizontal and vertical field of view (FOV) of a camera
%   given its focal length, pixel dimensions, and sensor resolution.
%
% INPUTS       Size     Type    Description                  Units
%   f          (1,1)    double  Focal length of the camera   [DU]
%   cam_W      (1,1)    double  Width resolution in pixels   [pixels]
%   cam_H      (1,1)    double  Height resolution in pixels  [pixels]
%   pix_W      (1,1)    double  Width of each pixel          [DU]
%   pix_H      (1,1)    double  Height of each pixel         [DU]
%
% OUTPUTS      Size     Type    Description                Units
%   FOV_W      (1,1)    double  Horizontal field of view   [rad]
%   FOV_H      (1,1)    double  Vertical field of view     [rad]
%
% NOTES
%
    
    % Converting width and height to physical dimensions
    W = cam_W * pix_W; % Camera width  [DU]
    H = cam_H * pix_H; % Camera height [DU]
    
    % Calculating field of view
    FOV_W = 2 * atan(W / (2 * f)); % Field of View Width  [rad]
    FOV_H = 2 * atan(H / (2 * f)); % Field of View Height [rad]
end