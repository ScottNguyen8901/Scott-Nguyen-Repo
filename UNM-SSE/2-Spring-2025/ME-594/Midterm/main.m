clear
close all
clc

%% =========================
%         USER INPUTS
% =========================

% ----- Observer / pass-finding -----
lat = 35.0844; lon = -106.6504; alt = 0;                 % [deg deg m]
tleFile    = 'ISS.tle';
startTime  = datetime(2022,11,17,0,0,0,'TimeZone','UTC');
endTime    = datetime(2022,11,18,0,0,0,'TimeZone','UTC');

sampleTime = 60;                                         % coarse step [s]
satElMin   = 10;                                         % [deg]
sunElMax   = -6;                                         % [deg]
padSec     = 30;                                         % refine padding [s]
fineStep   = 1;                                          % refine step [s]

% Choose which pass to simulate imagery for (after pass detection)
kPass = 1;                                               % <--- SELECT PASS INDEX HERE

% ----- Tracking error model (deterministic drift + jitter) -----
jitterAmp_deg = 0.0002;                                  % [deg] peak amplitude
jitterPerAz_s = 5;                                       % [s]
jitterPerEl_s = 6;                                       % [s]
driftRateAz_arcsec = 10;                                 % [arcsec/s]
driftRateEl_arcsec = 10;                                 % [arcsec/s]

% ----- Sensor / FOV (image sim) -----
Nx   = 1024;
Ny   = 1024;
FOVx = deg2rad(0.5);                                     % [rad]
FOVy = deg2rad(0.5);                                     % [rad]

% Time step for image generation (can be finer than fineStep)
dtImg = 0.1;                                             % [s] image integration step

% ----- PSF (Gaussian) -----
FWHM_pix = 2.0;                                          % [pix]

% ----- Photometry / noise params -----
% Intrinsic magnitude evolution (still allowed to vary):
m_obj_start = 12.0;                                      % mag at start of image timeline
m_obj_end   = 14.0;                                      % mag at end of image timeline (bigger = dimmer)

m_bg  = 21.0;                                            % sky brightness [mag/arcsec^2]

Phi0 = 1e10;                                             % photons/s/m^2 @ mag=0
A_ap = 0.05;                                             % m^2
eta  = 0.6;                                              % throughput (0-1)
readNoise_e = 5.0;                                       % e- rms

% ----- Playback -----
doPlayback  = true;
playbackFPS = 30;

% ----- Memory guard (optional) -----
storeFrames = true;                                      % if false, still generates but doesn't keep full cube

% ----- Noise realization control -----
nNoiseRealizations = 10;                                 % ALWAYS generate exactly 10 noise realizations

%% =========================
%   COARSE GRID (find passes)
% =========================

lla0 = [lat lon alt];
tle  = tleread(tleFile);

tC   = (startTime:seconds(sampleTime):endTime).';
utcC = [year(tC) month(tC) day(tC) hour(tC) minute(tC) second(tC)];
jdC  = juliandate(tC);
llaC = repmat(lla0, numel(tC), 1);

[rC_km, ~] = propagateOrbit(tC, tle);
rC_m   = (1e3*rC_km).';                                   % 3xN (m)
sunC_m = 1e3*planetEphemeris(jdC,'Earth','Sun','405','km');% Nx3 (m)

aerC_sat = eci2aer(rC_m,   utcC, llaC);                   % [az el rng] (deg,deg,m)
aerC_sun = eci2aer(sunC_m, utcC, llaC);

validC = (aerC_sat(:,2) >= satElMin) & (aerC_sun(:,2) <= sunElMax);
d      = diff([false; validC; false]);
iStart = find(d == 1);
iEnd   = find(d == -1) - 1;

fprintf('Found %d coarse pass window(s).\n', numel(iStart));

%% ===========================================
%   REFINE EACH PASS, ADD TRACK ERR + PHASE
% ===========================================

passes = cell(numel(iStart),1);
arcsec2deg = 1/3600;

for k = 1:numel(iStart)

    t0 = tC(iStart(k)) - seconds(padSec);
    tf = tC(iEnd(k))   + seconds(padSec);

    tF   = (t0:seconds(fineStep):tf).';
    utcF = [year(tF) month(tF) day(tF) hour(tF) minute(tF) second(tF)];
    jdF  = juliandate(tF);
    llaF = repmat(lla0, numel(tF), 1);

    [rF_km, ~] = propagateOrbit(tF, tle);
    rSat_m = 1e3*rF_km;                                   % Nx3 (m)
    rSat_m = rSat_m.';
    rSun_m = 1e3*planetEphemeris(jdF,'Earth','Sun','405','km'); % Nx3 (m)

    aerF_sat = eci2aer(rSat_m,   utcF, llaF);              % Nx3
    aerF_sun = eci2aer(rSun_m,   utcF, llaF);              % Nx3

    validF = (aerF_sat(:,2) >= satElMin) & (aerF_sun(:,2) <= sunElMax);
    idxV   = find(validF);

    if isempty(idxV)
        passes{k} = struct('tUTC',[],'jd',[],'aer_sat_true',[],'aer_sat_meas',[],'aer_sun',[], ...
                           'phaseDeg',[],'startUTC',NaT,'endUTC',NaT,'duration',0,'maxElSat',NaN, ...
                           'errAz_deg',[],'errEl_deg',[]);
        continue
    end

    tV  = tF(idxV);
    jdV = jdF(idxV);

    aer_true = aerF_sat(idxV,:);                           % [az el rng] deg,deg,m

    % --- deterministic jitter + linear drift (time since first valid sample)
    tau = seconds(tV - tV(1));                             % seconds from 0

    jitterAz_deg = jitterAmp_deg * sin(2*pi*tau/jitterPerAz_s);
    jitterEl_deg = jitterAmp_deg * cos(2*pi*tau/jitterPerEl_s);

    driftAz_deg  = (driftRateAz_arcsec * arcsec2deg) .* tau;
    driftEl_deg  = (driftRateEl_arcsec * arcsec2deg) .* tau;

    errAz_deg = jitterAz_deg + driftAz_deg;
    errEl_deg = jitterEl_deg + driftEl_deg;

    aer_meas = aer_true;
    aer_meas(:,1) = mod(aer_true(:,1) + errAz_deg, 360);
    aer_meas(:,2) = aer_true(:,2) + errEl_deg;

    % --- Solar phase angle (geocenter observer approximation)
    rSatV = rSat_m(idxV,:);                                % Nx3
    rSunV = rSun_m(idxV,:);                                % Nx3

    u_s2sun = rSunV - rSatV;
    u_s2obs = -rSatV;

    u_s2sun = u_s2sun ./ vecnorm(u_s2sun,2,2);
    u_s2obs = u_s2obs ./ vecnorm(u_s2obs,2,2);

    cang = sum(u_s2sun .* u_s2obs, 2);
    cang = max(min(cang,1),-1);
    phaseDeg = acosd(cang);

    passes{k} = struct( ...
        'tUTC',        tV, ...
        'jd',          jdV, ...
        'aer_sat_true',aer_true, ...
        'aer_sat_meas',aer_meas, ...
        'aer_sun',     aerF_sun(idxV,:), ...
        'phaseDeg',    phaseDeg, ...
        'startUTC',    tV(1), ...
        'endUTC',      tV(end), ...
        'duration',    seconds(tV(end) - tV(1)), ...
        'maxElSat',    max(aer_true(:,2)), ...
        'errAz_deg',   errAz_deg, ...
        'errEl_deg',   errEl_deg );
end

% --- basic pass summary
for k = 1:numel(passes)
    if isempty(passes{k}.tUTC), continue; end
    fprintf('Pass %d: %s -> %s | dur=%.1fs | maxEl=%.1f deg\n', ...
        k, string(passes{k}.startUTC), string(passes{k}.endUTC), ...
        passes{k}.duration, passes{k}.maxElSat);
end

% guard: pick a valid pass
if kPass < 1 || kPass > numel(passes) || isempty(passes{kPass}.tUTC)
    error('kPass=%d is invalid or empty. Choose a pass that has valid samples.', kPass);
end

%% ===========================================
%   IMAGE GENERATION for selected pass (kPass)
% ===========================================

P = passes{kPass};

% Sensor derived quantities
pixAngX = FOVx / Nx;                                      % [rad/pix]
pixAngY = FOVy / Ny;                                      % [rad/pix]

% Build image-time vector over the valid-window of the selected pass
t0 = P.tUTC(1);
tf = P.tUTC(end);
tImgUTC = (t0:seconds(dtImg):tf).';
tauImg  = seconds(tImgUTC - t0);                          % [s]

% Interpolate tracking errors onto image time (errors drive pixel offsets)
tauP = seconds(P.tUTC - t0);                              % pass time axis [s]
errAz_rad = deg2rad(interp1(tauP, P.errAz_deg, tauImg, 'spline', 'extrap'));
errEl_rad = deg2rad(interp1(tauP, P.errEl_deg, tauImg, 'spline', 'extrap'));

% Pixel location relative to center (tracker commands predicted center; residual = error)
xpix =  errAz_rad./pixAngX + Nx/2;
ypix = -errEl_rad./pixAngY + Ny/2;

%% ===========================================
%   Truncate to times where object is within FOV
%   - Supports MULTIPLE in-FOV segments (keeps longest)
%   - Builds fully-truncated struct P_FOV (top-level + sim + iod)
% ===========================================

inFOV = (xpix>=1 & xpix<=Nx & ypix>=1 & ypix<=Ny);
if ~any(inFOV)
    error('Object is never within the FOV for the selected pass.');
end

% Find ALL contiguous in-FOV segments on the image timeline
dF = diff([false; inFOV(:); false]);
segStart = find(dF==1);
segEnd   = find(dF==-1) - 1;
nSeg = numel(segStart);

% Keep the LONGEST segment
segLens = segEnd - segStart + 1;
[~, iKeep] = max(segLens);

idxIn  = segStart(iKeep);
idxOut = segEnd(iKeep);

tInUTC  = tImgUTC(idxIn);
tOutUTC = tImgUTC(idxOut);

% Map the image-time FOV window onto the PASS sample timeline (fineStep)
idxPassFOV = find(P.tUTC >= tInUTC & P.tUTC <= tOutUTC);
if isempty(idxPassFOV)
    error('FOV window did not overlap pass sample times. Check dtImg vs fineStep alignment.');
end

% Build fully-truncated pass struct
P_FOV = P;
P_FOV.tUTC         = P.tUTC(idxPassFOV);
P_FOV.jd           = P.jd(idxPassFOV);
P_FOV.aer_sat_true = P.aer_sat_true(idxPassFOV,:);
P_FOV.aer_sat_meas = P.aer_sat_meas(idxPassFOV,:);
P_FOV.aer_sun      = P.aer_sun(idxPassFOV,:);
P_FOV.phaseDeg     = P.phaseDeg(idxPassFOV);
P_FOV.errAz_deg    = P.errAz_deg(idxPassFOV);
P_FOV.errEl_deg    = P.errEl_deg(idxPassFOV);

P_FOV.startUTC  = P_FOV.tUTC(1);
P_FOV.endUTC    = P_FOV.tUTC(end);
P_FOV.duration  = seconds(P_FOV.endUTC - P_FOV.startUTC);
P_FOV.maxElSat  = max(P_FOV.aer_sat_true(:,2));

% Truncate IMAGE timeline arrays
tImgUTC   = tImgUTC(idxIn:idxOut);
tauImg    = tauImg(idxIn:idxOut);
xpix      = xpix(idxIn:idxOut);
ypix      = ypix(idxIn:idxOut);
errAz_rad = errAz_rad(idxIn:idxOut);
errEl_rad = errEl_rad(idxIn:idxOut);

Nimg  = numel(tImgUTC);
kStop = Nimg;

% Store sim fields
P_FOV.sim = struct();
P_FOV.sim.segments_img = [segStart(:), segEnd(:)];
P_FOV.sim.segment_kept = iKeep;
P_FOV.sim.idxIn_img    = idxIn;
P_FOV.sim.idxOut_img   = idxOut;
P_FOV.sim.kStop_img    = kStop;
P_FOV.sim.tImgUTC      = tImgUTC;
P_FOV.sim.tau_s        = tauImg;
P_FOV.sim.xpix         = xpix;
P_FOV.sim.ypix         = ypix;
P_FOV.sim.errAz_rad    = errAz_rad;
P_FOV.sim.errEl_rad    = errEl_rad;
P_FOV.sim.tInUTC       = tInUTC;
P_FOV.sim.tOutUTC      = tOutUTC;
P_FOV.sim.idxPassFOV   = idxPassFOV;

%% ===========================================
%   Time-varying angular rate -> time-varying magnitude
%   (motion-blur penalty model using omega0 derived from PSF + pixel scale)
% ===========================================

% Angular rates (rad/s) from pointing-error motion on the IMAGE timeline
dAz_dt = gradient(errAz_rad, dtImg);    % rad/s
dEl_dt = gradient(errEl_rad, dtImg);    % rad/s
omega_radps = hypot(dAz_dt, dEl_dt);    % rad/s

% omega0: rate where smear length over dtImg equals PSF FWHM (in pixels)
pixAng_min = min(pixAngX, pixAngY);     % rad/pix
omega0 = (FWHM_pix * pixAng_min) / dtImg;  % rad/s

% Motion-blur magnitude penalty:
delta_m = 2.5 * log10( sqrt(1 + (omega_radps./omega0).^2) );  % mag

% Intrinsic magnitude evolution (linear fade across the truncated image timeline)
tauEnd = max(tauImg);                                     % [s]
if tauEnd <= 0
    alpha = zeros(Nimg,1);
else
    alpha = tauImg ./ tauEnd;
end
m_intrinsic = m_obj_start + (m_obj_end - m_obj_start).*alpha;

% Final time-varying magnitude including smear penalty
m_obj_vec = m_intrinsic + delta_m;

% Save for later inspection
P_FOV.sim.omega_radps   = omega_radps;
P_FOV.sim.omega0_radps  = omega0;
P_FOV.sim.delta_m       = delta_m;
P_FOV.sim.m_intrinsic   = m_intrinsic;
P_FOV.sim.m_obj_vec     = m_obj_vec;

%% ===========================================
%   PSF (Gaussian)
% ===========================================

sigma = FWHM_pix/(2*sqrt(2*log(2)));
psfHalf = ceil(6*sigma);

[xg,yg] = meshgrid(-psfHalf:psfHalf);
PSF = exp(-(xg.^2 + yg.^2)/(2*sigma^2));
PSF = PSF/sum(PSF(:));

%% ===========================================
%   Photometry rates (time-varying)
% ===========================================

photonRate_obj_vec = Phi0*A_ap*eta*10.^(-0.4*m_obj_vec);        % photons/s

pixAngX_arcsec = pixAngX*180/pi*3600;
pixAngY_arcsec = pixAngY*180/pi*3600;
Omega_pix = pixAngX_arcsec*pixAngY_arcsec;                      % arcsec^2/pix
photonRate_bg_perPix = Phi0*A_ap*eta*10^(-0.4*m_bg)*Omega_pix;

%% ===========================================
%   Storage + noise (EXACTLY 10 noise realizations, but trajectory updates every dtImg)
% ===========================================

imgMean = zeros(Ny,Nx);

if storeFrames
    frames = zeros(Ny,Nx,kStop,'single');
else
    frames = [];
end

% Indices at which new noise is generated (exactly 10 times)
noiseIdx = unique(round(linspace(1, kStop, nNoiseRealizations)));
noiseIdx(1) = 1;

% Noise field that will be held constant between noiseIdx updates
noiseField = zeros(Ny,Nx);

for k = 1:kStop

    xc = round(xpix(k));
    yc = round(ypix(k));

    % Object PSF stamp into mean
    xIdx = (xc-psfHalf):(xc+psfHalf);
    yIdx = (yc-psfHalf):(yc+psfHalf);

    vX = xIdx>=1 & xIdx<=Nx;
    vY = yIdx>=1 & yIdx<=Ny;

    % Time-varying object brightness
    objRate = photonRate_obj_vec(k);                      % photons/s

    % Integrate object photons into mean
    imgMean(yIdx(vY),xIdx(vX)) = ...
        imgMean(yIdx(vY),xIdx(vX)) + PSF(vY,vX) * objRate * dtImg;

    % Integrate background into mean
    imgMean = imgMean + photonRate_bg_perPix * dtImg;

    % Update the NOISE FIELD only 10 times
    if any(k == noiseIdx)
        % Poisson approx: generate a noisy draw at current mean, then subtract mean -> pure noise field
        % (This preserves correct variance at this operating point.)
        poissonDraw = poissrnd(imgMean);
        readDraw    = readNoise_e * randn(Ny,Nx);

        noiseField = (poissonDraw - imgMean) + readDraw;
    end

    % Build current frame: mean + held noise field (trajectory updates every frame)
    imgFrame = imgMean + noiseField;
    imgFrame(imgFrame < 0) = 0;

    if storeFrames
        frames(:,:,k) = single(imgFrame);
    end
end

fprintf('Done generating %d frames.\n', kStop);

%% ===========================================
%   Playback AFTER generation (no zoom)
% ===========================================

if doPlayback && storeFrames
    figure('Color','w');
    hIm = imagesc(frames(:,:,1));
    axis image; set(gca,'YDir','normal');
    colormap hot; colorbar;
    xlabel('Pixel X'); ylabel('Pixel Y');

    mx = max(frames(:,:,end),[],'all');
    clim([0 max(50,0.9*mx)]);

    pausePlay = 1/playbackFPS;

    for k = 1:kStop
        set(hIm,'CData',frames(:,:,k));
        title(sprintf('Pass %d Playback | t = %.2f s | m = %.2f', kPass, tauImg(k), m_obj_vec(k)));
        drawnow limitrate;
        pause(pausePlay);
    end
elseif doPlayback && ~storeFrames
    fprintf('Playback skipped because storeFrames=false.\n');
end
