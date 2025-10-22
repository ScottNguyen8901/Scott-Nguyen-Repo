function [new_p, opt_pos_tau, opt_neg_tau] = Compute_p_for_bdry (bi_search, g, p_j, d_j, B, Delta)
%
% [new_p, opt_pos_tau, opt_neg_tau] = Compute_p_for_bdry (bi_search, g, p_j, d_j, B, Delta)
%  Searches for the minimum of the quadratuce function:
%      (g^t)*p + 0.5*(p^t)*B*p
%  where p = p_j + tau*d_j and ||p||=Delta.
%
% Input:
%   bi_search:      if (bi_search==1), both positive and negative values
%                   of tau are considered in the search.
%   g, p_j, d_j, B: define the quadratic (see above).
%
% Output:
%   new_p:                    holds the value of p.
%   opt_pos_tau, opt_neg_tau: hold the computed values for tau.
%
% Example:
%   bi_search = 1;
%   g   = [1; 2];
%   B   = [4  5; 7 8];
%   d_j = [4;  5];
%   p_j = [2; -1];
%   Delta = 5;
%   [new_p, opt_pos_tau, opt_neg_tau] = Compute_p_for_bdry (bi_search, g, p_j, d_j, B, Delta);
%   vec_norm(p_j + opt_pos_tau*d_j)- Delta
%   vec_norm(p_j + opt_neg_tau*d_j)- Delta


% For the Bisection method to work, we need:
% 1. Form the function for which we need the zero.
% 2. Pick two tau points: one for inside and one
%    for outside the trust region.


% Define function for computing ||p|| == Delta.
tau_fun = @(tau) vec_norm(p_j + tau*d_j) - Delta;

% Note that for tau=0, we have p_j that is always inside
% the trust region. We need a positive tau that takes us outside.
% This is accomplished using:
Small_tau = 0;
Large_tau = 1.2*Delta / abs(vec_norm(p_j) - vec_norm(d_j));
Max_iters = 100;

tols = [1.0e-10 1.0e-10];
[tau_vals, fun_evals] = bisection (Small_tau, Large_tau, tols, tau_fun, Max_iters);
opt_pos_tau = tau_vals(length(tau_vals));

if (bi_search == 0)
    new_p = p_j + opt_pos_tau*d_j;
    opt_neg_tau = 0;
    return;
end

% If we are doing bidirectional search, we must look for negative tau.
%  We can use the same Large_tau, but negative.
[tau_vals, fun_evals] = bisection (Small_tau, -Large_tau, tols, tau_fun, Max_iters);
opt_neg_tau = tau_vals(length(tau_vals));

% Evaluate between the two directions to see which one produces
% the minimum result.
p1 = p_j + opt_pos_tau*d_j;
p2 = p_j + opt_neg_tau*d_j;

if (((g.')*p1 + 0.5*(p1.')*B*p1) < (g.')*p2 + 0.5*(p2.')*B*p2)
  new_p = p1;
else
  new_p = p2;
end    

