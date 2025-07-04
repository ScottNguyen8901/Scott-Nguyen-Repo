function D = kepler_parabolic(d, p, dt, mu)
    MAXITER = 100; % in case the loop does not converge
    TOL = 1e-9;   % tolerance for convergence
    Dn = d;        % initial guess
    
    fn = (p*Dn + Dn^3/3) - 2*sqrt(mu)*dt;
    
    iter = 1; % count the iterations
    while (abs(fn) > TOL) && (iter <= MAXITER)
        Dn = Dn - fn / (p + Dn^2);
        fn = (p*Dn + Dn^3/3) - 2*sqrt(mu)*dt;
        iter = iter + 1;
        if iter >= MAXITER
            error('Did not converge');
        end
    end
    
    D = Dn; % converged solution
end
