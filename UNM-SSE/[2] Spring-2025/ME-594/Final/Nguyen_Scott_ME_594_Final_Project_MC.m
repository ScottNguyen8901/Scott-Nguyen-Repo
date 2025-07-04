clear
close all
clc

constants;

%% ------------------------Problem 2---------------------------------------

rng('shuffle');

% Initial setup
date0 = datetime('2025-05-05 03:18:34');
L   = deg2rad(35.0844);        % Albuquerque latitude
lon = deg2rad(-106.6504);    % Albuquerque longitude

% Initial ISS state in ECI [km; km/s]
r0 = [5165.8127; -1938.2678; -3951.6672];
v0 = [4.7525; 4.4652; 4.0249];
x_ISS_ECI_0 = [r0; v0];

% Time span and integration
tspan = 0:60*21;
opts = odeset('RelTol',1e-12, 'AbsTol',1e-12);
[t, x_ISS_ECI] = ode45(@(t,y) two_body_ode(t, y), tspan, x_ISS_ECI_0, opts);
r_ISS_ECI = x_ISS_ECI(:,1:3);

% Allocate arrays
N = length(t);
RA_Dec_true = zeros(N, 2); 
RA_Dec_pert = zeros(N, 2);
Az_El_true  = zeros(N, 2); 
Az_El_pert  = zeros(N, 2);

date_vec = date0 + seconds(t); % Date in seconds
sigma    = deg2rad(0.01);      % Noise std [rad]

num_runs    = 1000; % Number of MC runs
cond_thresh = 1e11; % Condition number threshold
cnd_cnt_itr = 0;    % Condition counter for each iteration
cnd_cnt_tot = 0;    % Condition counter for each MC run

% Allocate arrays
IOD_errors_all = zeros(num_runs, 6);
POD_errors_all = zeros(num_runs, 6);

for i = 1:num_runs
    
    clc
    percent_complete = (i / num_runs) * 100;
    fprintf('Run %d of %d complete (%.2f%%).\n', i, num_runs, percent_complete);
    for j = 1:N
        theta_g = GMST(date_vec(j));
        theta = lon + theta_g;
        Q = rot_ijk_sez(L, theta);
    
        r_th = Q * r_ISS_ECI(j,:)';
        u_th = (r_th - r_s_TH) / norm(r_th - r_s_TH);
    
        % Azimuth and Elevation
        az = mod(atan2(u_th(1), u_th(2)), 2*pi);
        el = asin(-u_th(3));
        Az_El_true(j,:) = [az, el];
        Az_El_pert(j,:) = [az, el] + sigma * randn(1,2);
    
        % RA and Dec
        r_s_eci = Q' * r_s_TH;
        u_eci = (r_ISS_ECI(j,:)' - r_s_eci) / norm(r_ISS_ECI(j,:)' - r_s_eci);
        RA = mod(atan2(u_eci(2), u_eci(1)), 2*pi);
        Dec = asin(u_eci(3));
        RA_Dec_true(j,:) = [RA, Dec];
        RA_Dec_pert(j,:) = [RA, Dec] + sigma * randn(1,2);
    end
    
    %------------------ Printing Tables and Performing IOD/POD ------------
    
    % Extract sample epochs and indices at t = 0, 60, 120 seconds (IOD) and every 60 seconds from t = 60 (POD)
    idx_IOD = [1, 61, 121];
    idx_POD = 121:60:length(t);
    
    % Extract measurement times
    date_meas_IOD = date_vec(idx_IOD);
    date_meas_POD = date_vec(idx_POD);
    
    % Extract RA/Dec and Az/El (true and perturbed) in radians for IOD
    RA_Dec_rad_true_IOD = RA_Dec_true(idx_IOD, :);
    RA_Dec_rad_pert_IOD = RA_Dec_pert(idx_IOD, :);
    Az_El_rad_true_IOD  = Az_El_true(idx_IOD, :);
    Az_El_rad_pert_IOD  = Az_El_pert(idx_IOD, :);
    
    % Extract RA/Dec and Az/El (true and perturbed) in radians for POD
    RA_Dec_rad_true_POD = RA_Dec_true(idx_POD, :);
    RA_Dec_rad_pert_POD = RA_Dec_pert(idx_POD, :);
    Az_El_rad_true_POD  = Az_El_true(idx_POD, :);
    Az_El_rad_pert_POD  = Az_El_pert(idx_POD, :);
    
    % Initial orbit estimation using Laplace method (IOD)
    x_ISS_ECI_IOD =...
        laplace_orbit_fit(r_s_TH, L, lon, date_meas_IOD, RA_Dec_rad_pert_IOD);
    
    % Initialize for POD estimation
    X0_i = x_ISS_ECI_IOD;  % Initial guess
    t0   = idx_POD(1);     % Start time for POD
    
    % Gauss-Newton parameters
    max_iter = 10;
    epsilon  = 1e-5;
    prev_RMS = 1e3;
    iter     = 1;
    
    % Gauss-Newton Loop
    while iter < max_iter
        H = [];  % Jacobian
        r = [];  % Residuals
    
        for k = 2:length(date_meas_POD)
            tf_POD = idx_POD(k);
            date_POD = date_meas_POD(k);
    
            % Perturbed Az/El measurement
            Y = Az_El_rad_pert_POD(k, :)';
            Y_hat =...
                compute_AzEl(X0_i, t0, tf_POD, r_s_TH, L, lon, date_POD);
            r_i = Y - Y_hat;
    
            % Jacobian via finite differencing
            H_i = zeros(2, length(X0_i));
            for l = 1:length(X0_i)
                dx = 0.01 * X0_i(l);
                X_plus  = X0_i; X_plus(l)  = X0_i(l) + dx;
                X_minus = X0_i; X_minus(l) = X0_i(l) - dx;
    
                az_el_plus  =...
                    compute_AzEl(X_plus,  t0, tf_POD, r_s_TH, L, lon, date_POD);
                az_el_minus =...
                    compute_AzEl(X_minus, t0, tf_POD, r_s_TH, L, lon, date_POD);
    
                H_i(:, l) = (az_el_plus - az_el_minus) / (2 * dx);
            end
    
            H = [H; H_i];
            r = [r; r_i];
        end
        
        cond_num = cond(H' * H);

        if cond_num > cond_thresh
            cnd_cnt_itr = cnd_cnt_itr + 1;
        end

        % Gauss-Newton update
        dX = (H' * H) \ (H' * r);
        X0_i = X0_i + dX;
    
        RMS = sqrt(r' * r / (2 * length(date_meas_POD)));
        fprintf('Iteration %d, RMS Residual = %.6e\n', iter, RMS);
    
        RMS_change = abs(RMS - prev_RMS) / prev_RMS;
        if RMS < epsilon || RMS_change < epsilon
            fprintf('Converged at iteration %d.\n', iter);
            break;
        end
    
        prev_RMS = RMS;
        iter = iter + 1;
    end
    
    if cnd_cnt_itr > 0
        cnd_cnt_tot = cnd_cnt_tot + 1;
    end

    % ---------------------- State Vector Comparison ----------------------
    
    % True and estimated states
    x_ISS_ECI_2   = x_ISS_ECI(61,:)';   % True state
    x_ISS_ECI_POD = X0_i;               % POD estimate
    
    IOD_error = x_ISS_ECI_2 - x_ISS_ECI_IOD;
    POD_error = x_ISS_ECI_2 - x_ISS_ECI_POD;
    
    IOD_errors_all(i, :) = IOD_error';
    POD_errors_all(i, :) = POD_error';

end

%% ----------------------------Plotting------------------------------------

fprintf(' %.2f%% (%d/%d) MC runs exceeded condition number threshold\n', ...
        100 * cnd_cnt_tot / num_runs, cnd_cnt_tot, num_runs);

% Define formatting options
fs = 12;  % Font size
fontname = 'Times New Roman';
fontweight = 'bold';
fontcolor = 'k';

% Component names for state vector
components = {'r_x', 'r_y', 'r_z', 'v_x', 'v_y', 'v_z'};

% Create histogram plot
figure;
for k = 1:6
    % IOD and POD subplots
    for i = 1:2
        subplot(6,2,2*k-1+i-1)
        errors = (i == 1) * IOD_errors_all(:,k) + (i == 2) * POD_errors_all(:,k);
        histogram(errors);
        
        % Title, labels, and styling
        title([components{k}, ' Error (', (i == 1) * 'IOD' + (i == 2) * 'POD', ')'], ...
            'FontName', fontname, 'FontSize', fs, 'FontWeight', fontweight, 'Color', fontcolor);
        xlabel('Error [km or km/s]', ...
            'FontName', fontname, 'FontSize', fs, 'FontWeight', fontweight, 'Color', fontcolor);
        ylabel('PDF', ...
            'FontName', fontname, 'FontSize', fs, 'FontWeight', fontweight, 'Color', fontcolor);
        set(gca, 'FontName', fontname, 'FontSize', fs, 'FontWeight', fontweight, 'XColor', fontcolor, 'YColor', fontcolor)
    end
end

% Super title
sgtitle('State Vector Error Histograms (IOD vs POD)', ...
    'FontName', fontname, 'FontSize', fs+2, 'FontWeight', fontweight, 'Color', fontcolor);