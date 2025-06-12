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
    
    R_pb = rot_1(theta);
    
    r = [(a + d) / 2;
         0;
         0];
    
    % Inertia tensors
    I_box_c_b = m_rp / 12 * diag([d^2 + h^2; ...
                               d^2 + h^2; ...
                               2 * d^2]);
    
    I_pan_c_p = m_p / 12 * diag([a^2; ...
                             2 * a^2; ...
                             a^2]);
    
    I_pan_c_b = R_pb' * I_pan_c_p * R_pb;

    % Inertia of the payload component in the center of mass frame
    I_pan_o_b = I_pan_c_b - m_p * (r' * r * eye(3) - r * r');
    I_box_o_b = I_box_c_b;

    % Total inertia of the system (spacecraft)
    I_sc_c = I_box_o_b + 2 * I_pan_o_b;
    
end