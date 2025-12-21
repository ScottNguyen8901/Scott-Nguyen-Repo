%% Drift-out-of-FOV simulation + PSF integration + background + noise
% FULL FRAME (no zoom)
% Mean integrated every dt
% Noise realized once per second
% Store frames
% Playback after generation

clear; clc; close all;

% ---------------------------
% Sensor / FOV configuration
% ---------------------------
Nx   = 1024;
Ny   = 1024;
FOVx = deg2rad(0.5);
FOVy = deg2rad(0.5);

pixAngX = FOVx / Nx;
pixAngY = FOVy / Ny;

% ---------------------------
% Time configuration
% ---------------------------
T  = 45;
dt = 0.1;
t  = (0:dt:T).';
N  = numel(t);

noiseDt   = 2.0;   % noise cadence [s]
nextNoise = 0.0;

% ---------------------------
% Tracking error model
% ---------------------------
jitterAz = deg2rad(0.0002) * sin(2*pi*t/5);
jitterEl = deg2rad(0.0002) * cos(2*pi*t/6);

driftRateAz_arcsec = 10;
driftRateEl_arcsec = 20;
arcsec2rad = (pi/180)/3600;

driftAz = driftRateAz_arcsec * arcsec2rad .* t;
driftEl = driftRateEl_arcsec * arcsec2rad .* t;

errAz = jitterAz + driftAz;
errEl = jitterEl + driftEl;

% ---------------------------
% Pixel location
% ---------------------------
xpix =  errAz./pixAngX + Nx/2;
ypix = -errEl./pixAngY + Ny/2;

inFOV = (xpix>=1 & xpix<=Nx & ypix>=1 & ypix<=Ny);
idxExit = find(~inFOV,1,'first');

if isempty(idxExit)
    kStop = N;
    fprintf('Object stays in FOV for %.1f s\n',T);
else
    kStop = idxExit;
    fprintf('Object exits FOV at t = %.2f s\n',t(idxExit));
end

% ---------------------------
% PSF (Gaussian)
% ---------------------------
FWHM_pix = 2.0;
sigma = FWHM_pix/(2*sqrt(2*log(2)));

psfHalf = ceil(6*sigma);
psfSize = 2*psfHalf + 1;

[xg,yg] = meshgrid(-psfHalf:psfHalf);
PSF = exp(-(xg.^2+yg.^2)/(2*sigma^2));
PSF = PSF/sum(PSF(:));

% ---------------------------
% Photometry
% ---------------------------
m_obj = 12.5;
m_bg  = 21.0;

Phi0 = 1e10;
A_ap = 0.05;
eta  = 0.6;

photonRate_obj = Phi0*A_ap*eta*10^(-0.4*m_obj);

pixAngX_arcsec = pixAngX*180/pi*3600;
pixAngY_arcsec = pixAngY*180/pi*3600;
Omega_pix = pixAngX_arcsec*pixAngY_arcsec;

photonRate_bg_perPix = Phi0*A_ap*eta*10^(-0.4*m_bg)*Omega_pix;

readNoise_e = 5.0;

% ---------------------------
% Storage
% ---------------------------
frames = zeros(Ny,Nx,kStop,'single');
imgMean = zeros(Ny,Nx);
imgNoisy_last = zeros(Ny,Nx);

fprintf('Generating frames: mean @ dt=%.2f s, noise @ %.1f s...\n',dt,noiseDt);

% ---------------------------
% Generate ALL frames (no plotting)
% ---------------------------
for k = 1:kStop

    if ~inFOV(k)
        kStop = k-1;
        frames = frames(:,:,1:kStop);
        break;
    end

    xc = round(xpix(k));
    yc = round(ypix(k));

    % Object PSF stamp into mean
    xIdx = (xc-psfHalf):(xc+psfHalf);
    yIdx = (yc-psfHalf):(yc+psfHalf);

    vX = xIdx>=1 & xIdx<=Nx;
    vY = yIdx>=1 & yIdx<=Ny;

    imgMean(yIdx(vY),xIdx(vX)) = ...
        imgMean(yIdx(vY),xIdx(vX)) + PSF(vY,vX)*photonRate_obj*dt;

    % Background integration
    imgMean = imgMean + photonRate_bg_perPix*dt;

    % Noise realization only once per second
    if t(k)>=nextNoise || k==1
        imgNoisy_last = poissrnd(imgMean);
        imgNoisy_last = imgNoisy_last + readNoise_e*randn(Ny,Nx);
        imgNoisy_last(imgNoisy_last<0)=0;
        nextNoise = nextNoise + noiseDt;
    end

    frames(:,:,k) = single(imgNoisy_last);
end

fprintf('Done generating %d frames.\n',kStop);

% ---------------------------
% Playback AFTER generation
% ---------------------------
playbackFPS = 30;                 % how fast the movie plays visually
pausePlay = 1/playbackFPS;

% figure('Color','w');
% hIm = imagesc(frames(:,:,1));
% axis image; set(gca,'YDir','normal');
% colormap hot; colorbar;
% xlabel('Pixel X'); ylabel('Pixel Y');
% 
% mx = max(frames(:,:,end),[],'all');
% clim([0 max(50,0.9*mx)]);

% for k = 1:kStop
%     set(hIm,'CData',frames(:,:,k));
%     title(sprintf('Playback | t = %.1f s', t(k)));
%     drawnow limitrate;
%     pause(pausePlay);
% end

%%
% ---------------------------
% Backtrack pixel offsets -> Az/El measurements
% ---------------------------

x_center = Nx/2;
y_center = Ny/2;

dx_pix = xpix(1:kStop) - x_center;
dy_pix = ypix(1:kStop) - y_center;

dAz =  dx_pix * pixAngX;     % rad
dEl = -dy_pix * pixAngY;     % rad

% Tracker believed pointing (example: perfect command)
Az_cmd = zeros(kStop,1);
El_cmd = zeros(kStop,1);

% Final measured angles
Az_meas = Az_cmd + dAz;
El_meas = El_cmd + dEl;
