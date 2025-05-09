function az_el = compute_AzEl(X0, t0, tf, r_s_TH, L, lon, date)
    [~, X_prop] = ode45(@two_body_ode, [t0 tf], X0);
    r_ISS = X_prop(end, 1:3)';

    [~, theta] = geocentric_to_ECI(r_s_TH, L, lon, date);
    Q_TH_ECI = rot_ijk_sez(L, theta);

    r_TH = Q_TH_ECI * r_ISS - r_s_TH;
    u_TH = r_TH / norm(r_TH);

    el = asin(-u_TH(3));
    az = mod(atan2(u_TH(1), u_TH(2)), 2*pi);

    if az < 0
        az = az + 2*pi;
    end

    az_el = [az; el];
end