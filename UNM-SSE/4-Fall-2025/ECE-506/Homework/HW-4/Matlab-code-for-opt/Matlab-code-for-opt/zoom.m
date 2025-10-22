function [alpha_opt] = zoom (phi_fun, phi_der_fun, phi0, phip0, ... 
                             a_lo, a_hi, phi_a_lo, phi_a_hi, phip_a_lo,  ...
                             c1, c2)
%
% [alpha_opt] = zoom (phi_fun, phi_der_fun, phi0, phip0, ... 
%                     a_lo, a_hi, phi_a_lo, phi_a_hi, phip_a_lo,  ...
%                     c1, c2)
% Computes the optimal value for alpha, so that the strong Wolfe conditions
% are satisfied. The procedure is an implmentation of the pseudocode given
% by Nocedal and Wright, page 60.
%
% Input:
%   phi_fun:     the phi-function to minimize
%   phi_der_fun: the derivative of phi to minimize 
%   c1, c2:      constants used in strong Wolfe conditions.
%   alpha_lo, alpha_hi: values for specifying alpha range. We must have:
%        1. Wolfe conditions satisfied between these two values.
%        2. phi(alpha_lo) < phi(alpha_hi) or other encoutered alphas.
%        3. phi'(alpha_lo)*(alpha_hi - alpha_lo) < 0.
%
% Output:
%   alpha_opt:   the value of alpha that satisfies the conditions.
%
% Example:
%  phi = @(x) 5*x^2 + 7*x - 2;
%  phi_der = @(x) 10*x + 7;
%  c1 = 1.0e-4;
%  c2 = 0.9;
%  phi0  = phi(0);
%  phip0 = phi_der(0);
%  a_lo  = -2;
%  a_hi  = +4;
%  phi_a_lo  = phi(a_lo);
%  phip_a_lo = phi_der(a_lo);
%  phi_a_hi  = phi(a_hi);
%  [alpha_opt] = zoom (phi, phi_der, phi0, phip0, ... 
%                      a_lo, a_hi, phi_a_lo, phi_a_hi, phip_a_lo,  ...
%                      c1, c2);
%  i=1;
%  for x=-5:0.1:5
%    phia(i)=phi(x);
%    i=i+1;
%  end
%  plot(-5:0.1:5, phia), hold on;
%  plot([0 a_lo a_hi], [phi0, phi_a_lo, phi_a_hi], '+');
%  plot(alpha_opt, phi(alpha_opt), '*');

max_iterations = 100;
k = 1;
while(k < max_iterations)
  if (abs(a_lo - a_hi)<1.0e-5)
      alpha_opt = a_lo;
      return;
  end
  
  % Use quadratic interpolation in the first iteration:  
  %  between a_lo and a_hi:  
  [a_j] = quad_step (a_lo, a_hi, phi_a_lo, phi_a_hi, phip_a_lo);
  
  phi_a_j = phi_fun(a_j); % phi(alpha_j)
  if ((phi_a_j > phi0 + c1*a_j*phip0) | (phi_a_j > phi_a_lo))
      % Generate new a_hi and phi(a_hi)
      a_hi = a_j;
      phi_a_hi = phi_a_j; 
  else
      % Generate new a_lo and phi'(a_lo).
      phip_a_j = phi_der_fun(a_j); % phi'(alpha_j)
      if (abs(phip_a_j) <= -c2*phip0)
          alpha_opt = a_j;
          return;
      end
      
      if (phip_a_j*(a_hi-a_lo) >= 0)
        a_hi = a_lo;
        
        % Update the function at a_hi:
        phi_a_hi = phi_a_lo;
      end    
      a_lo = a_j;
      
      % Update phi information at a_lo:
      phi_a_lo  = phi_a_j;
      phip_a_lo = phip_a_j;
  end    
  
  k = k+1;
end    

% If you are here, use a_lo!
alpha_opt = a_lo;