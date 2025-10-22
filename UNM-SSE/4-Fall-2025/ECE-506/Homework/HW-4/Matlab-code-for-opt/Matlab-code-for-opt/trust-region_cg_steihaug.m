function [p] = CG_Steihaug (vec_eps, g, B, Delta)
%  [p] = CG_Steihaug (vec_eps, g, B, Delta)
%
% Uses Conjugate-Gradient Steihaug's approach to compute
% a direction p that minimizes the quadratic function defined 
% by B and initial_r. The code is taken from Nocedall and Wright,
% page 75.
%
% Input:
%   vec_eps:   terminate if residual < vec_eps*||r0||
%   initial_r: initial value for the residual.
%   B:         an estimate of the Hessian 
% Output:
%   p: the direction minimizing the quadratic.
%
% Example:
% vec_eps = 1.0e-4;
% fun   = @(x) x(1).^2 + 5*x(2).^2;
% g_fun = @(x) [2*x(1); 10*x(2)];
% Hessian_fun = @(x) [2, 0; 0 10];
% initial_x = [-2; -3];
% g = g_fun(initial_x);
% B = Hessian_fun(initial_x);
% Delta = 100;
% [p] = CG_Steihaug (vec_eps, g, B, Delta);
% disp('Function at x+p');
% fun(initial_x + p)
% disp('Gradient at x+p');
% g_fun(initial_x + p)

% Initialize the j-th variables for j=0:
r_j = g;
r0  = g;
p_j = zeros(length(g), 1);
d_j = -r_j;

% Test for r0=0:
if (vec_norm(r_j)<vec_eps)
  p = p_j;
  return;
end

% Compute all directions
for j=0:1:length(g)-1
    % Negative curvature direction?
    if ((d_j .')*B*d_j <= 0)
        disp('Negative direction found.');
        bi_search = 1;
        [new_p, opt_pos_tau, opt_neg_tau] = Compute_p_for_bdry (bi_search, g, p_j, d_j, B, Delta);
        p = new_p;
        return;
    end        
    
    a_j        = ((r_j .')*r_j)/((d_j .')*B*d_j);
    p_j_plus_1 = p_j + a_j*d_j;
    
    % If we are beyond the trust region, reduce vector length:
    if (vec_norm(p_j_plus_1)>Delta)
        disp('Reached the trust region boundary');
        bi_search = 0; % Only search for positive tau.
        [new_p, opt_pos_tau, opt_neg_tau] = Compute_p_for_bdry (bi_search, g, p_j, d_j, B, Delta);
        p = new_p;
        return;
    end        
        
    r_j_plus_1 = r_j + a_j*B*d_j;
    if (vec_norm(r_j_plus_1) < vec_eps*vec_norm(r0))
        disp('Met Stopping Criterion');
        p = p_j_plus_1;
        return;
    end
        
    % Update:
    beta_j_plus_1 = ((r_j_plus_1 .') * r_j_plus_1) / ((r_j .') * r_j);
    d_j_plus_1    = -r_j_plus_1 + beta_j_plus_1*d_j;
    
    % Copy current (j+1)-th variable as j-th variable.
    r_j = r_j_plus_1;
    p_j = p_j_plus_1;
    d_j = d_j_plus_1;    
end    

% Check if you are here w/out a valid direction:
disp('Examined all directions!');