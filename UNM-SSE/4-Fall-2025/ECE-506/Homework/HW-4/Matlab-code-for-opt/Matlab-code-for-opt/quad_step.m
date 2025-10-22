function [a_j] = quad_step (a_lo, a_hi, phi_a_lo, phi_a_hi, phip_a_lo)
%
% [a_j] = quad_step (a_lo, a_hi, phi_a_lo, phi_a_hi, phip_a_lo)
%  fits a quadratic model to compute the point where phi is minimum,
%  subject to the condition that the estimated a is between a_lo and a_hi.
%
% Inputs:
%   a_lo, a_hi:          the alpha values determining the range for alpha.
%   phi_a_lo, phi_a_hi:  phi(a_lo), phi(a_hi)
%   phip_a_lo:           phi'(a_lo)
%
% Outputs:
%   a_j: the estimated value where the function is minimized.
%
% Example:
%  phi = @(x) 5*x^2 + 7*x - 2;
%  phi_der = @(x) 10*x + 7;
%  a_lo = -2;
%  phi_a_lo  = phi(a_lo);
%  phip_a_lo = phi_der(a_lo);
%  a_hi = 1;
%  phi_a_hi = phi(a_hi);
%  [a_j] = quad_step (a_lo, a_hi, phi_a_lo, phi_a_hi, phip_a_lo);
%  i=1;
%  for x=-5:0.1:5
%    phia(i)=phi(x);
%    i=i+1;
%  end
%  plot([a_lo a_hi], [phi_a_lo phi_a_hi], '+'), hold on;
%  plot(a_j, phi(a_j), '*');
%  plot(-5:0.1:5, phia);


% Apply quadratic interpolation using:
%  a_lo, phi_a_lo

% For quadratic interpolation, use:
%   phi(x) = a*x^2 + b*x + c.
% At x=a_lo, a_hi:
%   phi(a_lo) = a*(a_lo)^2 + b*a_lo + c (1)
%   phi(a_hi) = a*(a_hi)^2 + b*a_hi + c (2)
% For the derivative at x=a_lo:
%   phi'(a_lo) = 2*a*a_lo + b   (3)
% This allows us to setup:
%  [2*a_lo 1     0][a] = [phi'(a_lo)]
%  [a_lo^2 a_lo  1][b]   [phi (a_lo)]
%  [a_hi^2 a_hi  1][c]   [phi (a_hi)]
% For the minimum point, we require:
%  phi'(x) = 2*a*x + b = 0
%  => x = -b/(2*a)
% phi''(x) = 2*a
%

% Set up the system of equations:
A    = [2*a_lo 1 0; a_lo^2 a_lo 1; a_hi^2 a_hi 1];
Rhs  = [phip_a_lo; phi_a_lo; phi_a_hi];
Soln = A \ Rhs;

% Extract the quadratic coefficients:
a = Soln(1);
b = Soln(2);
c = Soln(3);

if (a<0)
  if (phi_a_lo < phi_a_hi)  
    a_j = a_lo; % Cannot do better!  
    return;
  else
    a_j = a_hi;
    return;
  end  
end    

a_j = -b/(2*a);
if (a_lo < a_hi)
  if ((a_j >= a_lo) & (a_j <= a_hi))
      return;
  else % Pick an endpoint:   
      if (phi_a_lo <= phi_a_hi)
        a_j = a_lo;
        return;
      else
        a_j = a_hi;
        return;
      end  
  end
else % a_lo >= a_hi
  if ((a_j >= a_hi) & (a_j <= a_lo))
      return;
  else
      if (phi_a_lo <= phi_a_hi)
        a_j = a_lo;
        return;
      else
        a_j = a_hi;
        return;
      end
  end
end


