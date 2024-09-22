function E = kepler(M, e)
    max_iter = 100; % in case the loop does not converge
    tol = 1e-11; % tolerance for convergence
    E_n = M; % first guess
    f_n = E_n - e*sin(E_n) - M; % f(E) that we want to be zero
    iter = 1; % count the iterations
    while ( abs(f_n) > tol ) & iter > tol
        E_n = E_n - f_n/(1-e*cos(E_n));
        f_n = E_n - e*sin(E_n) - M;
        iter = iter+1;
        if (iter>=max_iter)
        error('Did not converge');
        end
    end
    E = E_n;
end
