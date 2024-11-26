function I_sc_c = comp_intertia(theta, m_rp, d, h, a, rho)
    % DESCRIPTION
    %   Compute the inertia matrix of a spacecraft system 
    %
    % INPUTS        size    Type
    %   theta       (1,1)   Double  Angle of rotation [degrees]
    %   m_rp        (1,1)   Double  Mass of bus [kg]
    %   m_p         (1,1)   Double  Mass of solar panal[kg]
    %   d           (1,1)   Double  Base length of bus[m]
    %   h           (1,1)   Double  Height of bus [m]
    %   a           (1,1)   Double  Length/width of solar panal[m]
    %
    % OUTPUTS       size    Type
    %   I_sc_c      (3,3)   Double  Total inertia matrix of the spacecraft
    %
    % FUNCTION
    %
    
    A = a^2;   
    m_p = rho * A; 
    
    R_1 = rot_1(theta);
    
    r = [(d + a) / 2;
         0;
         0];
    
    % Inertia tensors
    I_rp_c = m_rp / 12 * diag([d^2 + h^2; ...
                                2 * d^2; ...
                                d^2 + h^2]);
    
    I_p_o = m_p / 12 * diag([a^2; ...
                             a^2; ...
                             2 * a^2]);
    
    % Inertia of the payload component in the center of mass frame
    I_p_c = R_1 * I_p_o * R_1' - m_p * (r' * r * eye(3) - r * r');
    
    % Total inertia of the system (spacecraft)
    I_sc_c = I_rp_c + 2 * I_p_c;
    
end
