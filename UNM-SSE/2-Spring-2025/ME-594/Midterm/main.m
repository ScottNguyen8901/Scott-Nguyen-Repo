clear
close all
clc

%% =========================
%         USER INPUTS
% =========================

% Observer
lat = 35.0844;
lon = -106.6504;
alt = 0;                              % [deg deg m]

tleFile   = 'ISS.tle';

startTime = datetime(2026,11,17,0,0,0,'TimeZone','UTC');
endTime   = datetime(2026,11,18,0,0,0,'TimeZone','UTC');

% Pass detection
sampleTime = 60;                      % coarse step [s]
fineStep   = 1;                       % refined step [s]
padSec     = 30;                      % padding [s]

satElMin = 10;                        % [deg]
sunElMax = -6;                        % [deg]

kPass = 1;                            % select pass index

%% =========================
%   STATION ECEF POSITION
% =========================

wgs84 = wgs84Ellipsoid('meter');
[xSta_ecef, ySta_ecef, zSta_ecef] = geodetic2ecef(wgs84, lat, lon, alt);
rSta_ecef_m = [xSta_ecef; ySta_ecef; zSta_ecef];   % 3x1 [m]

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
rC_m   = (1e3*rC_km).';                                    % 3xN (m)
sunC_m = 1e3*planetEphemeris(jdC,'Earth','Sun','405','km');

aerC_sat = eci2aer(rC_m,   utcC, llaC);                    % [az el rng]
aerC_sun = eci2aer(sunC_m, utcC, llaC);

validC = (aerC_sat(:,2) >= satElMin) & ...
         (aerC_sun(:,2) <= sunElMax);

d      = diff([false; validC; false]);
iStart = find(d == 1);
iEnd   = find(d == -1) - 1;

fprintf('Found %d coarse pass window(s).\n', numel(iStart));

%% =====================================
%   REFINE PASSES (AZ / EL / TIME ONLY)
%   + INCLUDE STATION ECEF IN STRUCT
% =====================================

passes = cell(numel(iStart),1);

for k = 1:numel(iStart)

    t0 = tC(iStart(k)) - seconds(padSec);
    tf = tC(iEnd(k))   + seconds(padSec);

    tF   = (t0:seconds(fineStep):tf).';
    utcF = [year(tF) month(tF) day(tF) hour(tF) minute(tF) second(tF)];
    jdF  = juliandate(tF);
    llaF = repmat(lla0, numel(tF), 1);

    [rF_km, ~] = propagateOrbit(tF, tle);
    rSat_m = (1e3*rF_km).';                                % 3xN (m)
    rSun_m = 1e3*planetEphemeris(jdF,'Earth','Sun','405','km');

    aerF_sat = eci2aer(rSat_m, utcF, llaF);
    aerF_sun = eci2aer(rSun_m, utcF, llaF);

    validF = (aerF_sat(:,2) >= satElMin) & ...
             (aerF_sun(:,2) <= sunElMax);

    idxV = find(validF);
    if isempty(idxV)
        passes{k} = struct([]);
        continue
    end

    tV  = tF(idxV);
    jdV = jdF(idxV);

    az_deg = aerF_sat(idxV,1);
    el_deg = aerF_sat(idxV,2);

    passes{k} = struct( ...
        'tUTC',        tV, ...
        'jd',          jdV, ...
        'az_deg',      az_deg, ...
        'el_deg',      el_deg, ...
        'rSta_ecef_m', rSta_ecef_m, ...
        'lat_deg',     lat, ...
        'lon_deg',     lon, ...
        'alt_m',       alt, ...
        'startUTC',    tV(1), ...
        'endUTC',      tV(end), ...
        'duration',    seconds(tV(end)-tV(1)), ...
        'maxEl_deg',   max(el_deg) );
end

%% =========================
%       PASS SUMMARY
% =========================

for k = 1:numel(passes)
    if isempty(passes{k}), continue; end
    fprintf('Pass %d: %s -> %s | dur=%.1fs | maxEl=%.1f deg\n', ...
        k, string(passes{k}.startUTC), string(passes{k}.endUTC), ...
        passes{k}.duration, passes{k}.maxEl_deg);
end

%% =========================
%       GUARD CHECK
% =========================

if kPass < 1 || kPass > numel(passes) || isempty(passes{kPass})
    error('kPass=%d is invalid or empty.', kPass);
end
