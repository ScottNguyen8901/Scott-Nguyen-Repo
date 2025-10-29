close all
clear all

iters = [5, 8, 12, 20];
for iter=1:length(iters)
    N = iters(iter);
    
    % Build the matrix
    [j,i] = meshgrid(1:N, 1:N);
    A  = 1 ./ (i+j-1);
    b  = ones(N,1);
    x0 = zeros(N,1);
    vec_eps = 1e-6;

   [x_opt, num_of_iters] = conj_grad (x0, A, b, vec_eps);

   disp(sprintf('N=%d, iterations=%d', N, num_of_iters));
   disp(sprintf('Residual norm = %10.8f < %10.8f', ...
        vec_norm(A*x_opt - b), vec_eps));
end    
