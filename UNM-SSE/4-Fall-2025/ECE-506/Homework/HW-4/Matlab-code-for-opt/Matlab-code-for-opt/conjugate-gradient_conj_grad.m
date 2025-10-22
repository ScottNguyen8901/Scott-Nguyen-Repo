function [x_opt, num_of_iters] = conj_grad (initial_x, A, b, vec_eps)
%
% [x_opt, num_of_iters] = conj_grad (initial_x, A, b, vec_eps)
% Solves Ax=b using the conjugate gradient method, up to
% a residual accuracy spevified by vec_eps.
%
% Input:
%  initial_x: starting value for x
%  A, b:      matrix A, vector b, so that Ax=b
%  vec_eps:   residual tolerance.
%
% Output:
%  x_opt:     the estimated value for x.
%
% Example:
%  A = [1, 2; 3 0];
%  b = [3; 4];
%  vec_eps   = 0.00001;
%  initial_x = [0; 0];
%  [x_opt, num_of_iters] = conj_grad (initial_x, A, b, vec_eps)
%  A*x_opt - b

% Set the values for k=0:
x_k = initial_x;
r_k = A*x_k - b;
p_k = -r_k;
k   = 0;

while (vec_norm(r_k) > vec_eps) 
  alpha_k    = ((r_k.')*r_k)/((p_k.')*A*p_k);
  x_k_plus_1 = x_k + alpha_k*p_k;
  
  r_k_plus_1 = r_k + alpha_k*A*p_k;
  
  beta_k_plus_1 = ((r_k_plus_1.')*(r_k_plus_1)) / ((r_k.')*r_k);
  p_k_plus_1    = -r_k_plus_1 + beta_k_plus_1*p_k;
  
  % Update the previous values:
  x_k = x_k_plus_1;
  r_k = r_k_plus_1;
  p_k = p_k_plus_1;
  
  k = k+1;
end    
x_opt = x_k;
num_of_iters = k;