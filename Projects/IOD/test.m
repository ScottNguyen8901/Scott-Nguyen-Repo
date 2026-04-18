clear; clc;

%% Time
startTime = datetime(2026,1,15,0,0,0,'TimeZone','UTC');
dt = 10;
duration = 2*3600;
time = startTime + seconds(0:dt:duration);
m = numel(time);

%% Orbit (meters + degrees) using propagateOrbit
a_m  = 7080.6e3;
ecc  = 0.0015;
incl = 98.20;
RAAN = 95.21;
argp = 120.48;
nu0  = 10.0;

pos_icrf_m = propagateOrbit(time, a_m, ecc, incl, RAAN, argp, nu0, ...
    PropModel="two-body-keplerian");    % 3-by-m
r_icrf = pos_icrf_m.';                  % m-by-3

%% New Mexico station (Albuquerque example)
gsLat_deg = 35.0844;
gsLon_deg = -106.6504;
gsAlt_m   = 1620;

gs_ecef_m = lla2ecef([gsLat_deg, gsLon_deg, gsAlt_m]);
gs_ecef_m = gs_ecef_m(:).';

%% WGS84 spheroid for Mapping Toolbox transforms
wgs84 = wgs84Ellipsoid("meter");

%% Measurements (no noise)
az_deg  = nan(m,1);
el_deg  = nan(m,1);
ra_deg  = nan(m,1);
dec_deg = nan(m,1);
u_icrf  = nan(m,3);
range_m = nan(m,1);

for k = 1:m
    % Station ECEF -> ICRF
    gs_icrf_m = ecef2eci(time(k), gs_ecef_m);
    gs_icrf_m = gs_icrf_m(:).';

    % LOS in ICRF
    rho_icrf = r_icrf(k,:) - gs_icrf_m;
    range_m(k) = norm(rho_icrf);
    u = rho_icrf / range_m(k);
    u_icrf(k,:) = u;

    % RA/Dec
    ra = atan2d(u(2), u(1)); if ra < 0, ra = ra + 360; end
    dec = asind(u(3));
    ra_deg(k)  = ra;
    dec_deg(k) = dec;

    % Satellite ICRF -> ECEF
    r_ecef_m = eci2ecef(time(k), r_icrf(k,:));
    r_ecef_m = r_ecef_m(:).';

    rho_ecef = r_ecef_m - gs_ecef_m;

    % Mapping Toolbox: ECEF -> NED (needs spheroid + degrees flag)
    [n, e, d] = ecef2ned( ...
        rho_ecef(1), rho_ecef(2), rho_ecef(3), ...
        gsLat_deg, gsLon_deg, gsAlt_m, ...
        wgs84, true);

    % Convert NED -> ENU then ENU -> AER using enu2aer
    xEast  = e;
    yNorth = n;
    zUp    = -d;

    [az, el, ~] = enu2aer(xEast, yNorth, zUp, "degrees");
    az_deg(k) = az;
    el_deg(k) = el;
end

visible = el_deg > 0;

meas.time    = time(:);
meas.az_deg  = az_deg;
meas.el_deg  = el_deg;
meas.ra_deg  = ra_deg;
meas.dec_deg = dec_deg;
meas.u_icrf  = u_icrf;
meas.range_m = range_m;
meas.visible = visible;

disp("Done. Struct 'meas' created.");
fprintf("Visible samples: %d / %d\n", nnz(visible), m);
