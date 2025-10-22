function [alpha_interp] = interp_cubic (alpha0, alpha1, phi0, phia0, phia1, phip0)
%
% [alpha_interp] = interp_cubic (alpha0, alpha1, phi0, phia0, phia1, phip0)
%  Compute alpha_interp where phi(alpha) is minimized, using cubic
%  interpolation.
%
% Input:
%  phi0, phip0:    original values of phi and phi' at the origin.
%  alpha0, alpha1: alpha values where phi is known.
%  phia0, phia1:   corresponding values for phi.
%
% Output:
%  alpha_interp:   an interpolated value for alpha.
%
% Example:
%  phi = @(x) 5*x^3 + 5*x^2 + 7*x - 2;
%  phi_der = @(x) 15*x^2 + 10*x + 7;
%  phi0   = phi(0);
%  phip0  = phi_der(0);
%  alpha0 = 1;
%  alpha1 = 2;
%  phia0  = phi(alpha0);
%  phia1  = phi(alpha1);
%  [alpha_interp] = interp_cubic (alpha0, alpha1, phi0, phia0, phia1, phip0);
%  i=1;
%  for x=-5:0.1:5
%    phia(i)=phi(x);
%    i=i+1;
%  end
%  plot([0 alpha0 alpha1], [phi0 phia0 phia1], '.'), hold on;
%  plot(alpha_interp, phi(alpha_interp), '*');
%  plot(-5:0.1:5, phia);

% Cubic fit:
const_mult = 1/(alpha0^2*alpha1^2*(alpha1-alpha0));
A_matrix   = [alpha0^2  -alpha1^2; ...
              -alpha0^3  alpha1^3];
b_vector   = [phia1 - phi0 - phip0*alpha1; ...
              phia0 - phi0 - phip0*alpha0];
a_b_vec = const_mult*A_matrix*b_vector;

a = a_b_vec(1);
b = a_b_vec(2);

root_val = b^2 - 3*a*phip0;
if (root_val >= 0)
  alpha2 = (-b + sqrt(root_val))/(3*a);
else
  alpha2 = -b/(3*a);  
end    
alpha_interp = alpha2;