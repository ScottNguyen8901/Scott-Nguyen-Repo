function [alpha_interp] = interp_quad (alpha0, phi0, phip0, phia0)
%
%  [alpha_interp] = interp_quad (alpha0, phi0, phip0, phia0)
%  Computes alpha where phi(a) is minimized, using quadratic interpolation.
%
% Inputs:
%   alpha_0: the new point where phi is known.
%   phia0:   the corresponding phi(alpha_0).
%   phip0:   original phi'(0)
%   phi0:    original phi(0)
% Outputs:
%   alpha_interp: is the value of alpha
%
% Example:
%  phi = @(x) 5*x^2 + 7*x - 2;
%  phi_der = @(x) 10*x + 7;
%  phi0   = phi(0);
%  phip0  = phi_der(0);
%  alpha0 = 1;
%  phia0  = phi(alpha0);
%  [alpha_interp] = interp_quad (alpha0, phi0, phip0, phia0);
%  i=1;
%  for x=-5:0.1:5
%    phia(i)=phi(x);
%    i=i+1;
%  end
%  plot([0 alpha0], [phi0 phia0], '.'), hold on;
%  plot(alpha_interp, phi(alpha_interp), '*');
%  plot(-5:0.1:5, phia);

% Quadratic fit:
alpha_interp = -phip0*alpha0^2 / (2*(phia0 - phi0 - phip0*alpha0));
