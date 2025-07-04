function rho = comp_rho(h, atmo_data)
    % Computes atmospheric density [kg/m^3] given altitude [km] and data.
    % 
    % INPUTS
    %   h           (1,1)    Double     Altitude above Earth's surface         [km]
    %   atmo_data   (N,4)    Double     Atmospheric data matrix containing:     []
    %                                   [altitude_max, base_alt, rho0, H]       [km, km, kg/m^3, km]
    %
    % OUTPUTS
    %   rho_km      (1,1)    Double     Atmospheric density [kg/m^3]           [kg/m^3]

    if h < 0
        rho = 0;
    else
        % Find the index corresponding to the closest altitude in the atmo_data
        idx = find(atmo_data(:, 1) >= h, 1, 'first');
        
        if isempty(idx)
            rho = 0;  % If no data for the altitude, set density to 0
        else
            h0 = atmo_data(idx, 2);   % Base altitude [km]
            rho0 = atmo_data(idx, 3); % Base density [kg/m^3]
            H = atmo_data(idx, 4);    % Scale height [km]
            rho = rho0 * exp(-(h - h0) / H);  % [kg/m^3]
        end
    end
end