function [p_dogleg] = dogleg (g, B, Delta)
% 
% [p_dogleg] = dogleg (g, B, Delta)
% Applies the Dogleg method (page 72 of Nocedal & Wright) to
% compute a vector p that minimizes the quadratic defined
% in terms of the gradient g and the Hessian B.
%
% Inputs:
%   g:     the value of the gradient at a set point.
%   B:     the value of the Hessian at a set point.
%   Delta: the trust-region length in radius.
%
% Outputs:
%   p_dogleg: the dogleg solution to the problem. 
%
% Example:
% fun   = @(x) x(1).^2 + 5*x(2).^2;
% g_fun = @(x) [2*x(1); 10*x(2)];
% Hessian_fun = @(x) [2, 0; 0 10];
% initial_x = [-2; -3];
% g = g_fun(initial_x);
% B = Hessian_fun(initial_x);
% Delta = 3.4;
% [p_dogleg] = dogleg (g, B, Delta);
% disp('Function at x+p');
% fun(initial_x + p_dogleg)
% disp('Gradient at x+p');
% g_fun(initial_x + p_dogleg)
%

p_U = -((g.')*g)/((g.')*B*g) * g;
p_B = -B\g; 

% For small Delta, pick steepest descent:
p_U_norm = vec_norm(p_U);
if (Delta <= p_U_norm)
   p_dogleg = (Delta/p_U_norm) * p_U; 
   return;
end    

% For large Delta, simply use the Newton minimum:
p_B_norm = vec_norm(p_B);
if (Delta >= p_B_norm)
   p_dogleg = p_B;
   return;
end    

% If we are here, then the vector norm is between
% the previous two. We must use Bisection to find the
% minimum.
bi_search = 0; % For the dogleg, we only consider positive values for tau.
p_j = p_U;     % Search for tau in  p_j + tau*d_j
d_j = p_B - p_U; 
Delta
vec_norm(p_U)
vec_norm(p_B)
[new_p, opt_pos_tau, opt_neg_tau] = Compute_p_for_bdry (bi_search, g, p_j, d_j, B, Delta);

opt_pos_tau
p_dogleg = new_p;
