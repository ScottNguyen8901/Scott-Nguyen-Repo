function F = kepler_hyperbolic(M, e)
    MAXITER = 100; % maximum number of iterations
    TOL = 1e-11; % tolerance for convergence
    Fn = (6*M)^(1/3); % initial guess
    fn = e * sinh(Fn) - Fn - M; % function value
    iter = 1; % iteration counter
    
    while abs(fn) > TOL
        % Update guess using Newton's method
        Fn = Fn - fn / (e * cosh(Fn) - 1);
        % Update function value
        fn = e * sinh(Fn) - Fn - M;
        % Increment iteration counter
        iter = iter + 1;
        % Check for maximum iterations
        if iter >= MAXITER
            error('Did not converge');
        end
    end
    
    % Return the solution
    F = Fn;
end
