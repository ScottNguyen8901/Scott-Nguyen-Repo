function [A, B, Beta] = gen_HESO(m, B0, U, w_o)
    % DESCRIPTION
    %   Computes the matrices A, B, and Beta for the High-Order Extended 
    %   State Observer (HESO) based on the given parameters for MIMO systems.
    %
    % INPUTS       size         Type        Description
    %   m          (1, 1)       Integer     Number of states for the observer
    %   B0         (m, p)       Double      Control input matrix (m states, p inputs)
    %   U          (p, 1)       Double      Control input vector (p inputs)
    %   w_o        (1, 1)       Double      Constant parameter for observer gains
    %
    % OUTPUTS      size         Type        Description
    %   A          (m*p, p)     Double      State matrix for the observer
    %   B          (m*p, p)     Double      Control influence matrix for the observer
    %   Beta       (m*p, m*p)   Double      Observer gain matrix (Beta)
    %
    % FUNCTION

    % Step 1: Generate matrix A (state matrix)
    A = zeros(m * p, p);  % Initialize A as an (m*p)-by-p matrix
    for i = 1:m*p - p
        A(i, i + p) = 1;  % Set the superdiagonal elements to 1 for coupling
    end

    % Step 2: Generate matrix B (control influence matrix)
    B = zeros(m * p, p);  % Initialize B as an (m*p)-by-p matrix
    for i = 1:m
        B((i-1)*p+1:i*p, :) = B0;  % Set blocks of B to B0
    end

    % Step 3: Generate matrix Beta (observer gain matrix)
    Beta = zeros(m * p, m * p);  % Initialize Beta as (m*p)-by-(m*p) matrix
    for i = 1:m * p
        % Using the formula for β_i with w_o^i
        Beta(i, i) = factorial(m * p) / (factorial(i) * factorial(m * p - i)) * w_o^(i-1);
    end
end