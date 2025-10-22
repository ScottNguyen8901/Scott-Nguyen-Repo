function [a_opt] = line_search (phi_fun, phi_der_fun, phi0, phip0, initial_a_i, a_max, c1, c2)
%
% [a_opt] = line_search (phi_fun, phi_der_fun, phi0, phip0, initial_a_i, a_max, c1, c2)
% Performs a line-search to compute an a_i that satisfies the strong
% Wolfe conditions. It implements Algorithm 3.2 page 59 of Nocedal &
% Wright.
%
% Inputs:
%   phi_fun:     a function for evaluating the local model.
%   phi_der_fun: a function for evaluating the derivative of the model.
%   phi0, phip0: initial function values phi(0),phi'(0).
%   initial_a_i: allows us to initialize the value for a_i, (See page 60).
%   a_max:       is the maximum value allowed for alpha.
%
% Outputs:
%   a_opt:   optimal value for alpha that satisfies the strong Wolfe
%                conditions.
% Example:
%  phi = @(x) 5*x^2 - 15*x - 2;
%  phi_der = @(x) 10*x - 15;
%  c1 = 1.0e-4;
%  c2 = 0.9;
%  phi0  = phi(0);
%  phip0 = phi_der(0);
%  initial_a_i = 5;
%  a_max       = 10;
%  [a_opt] = line_search (phi, phi_der, phi0, phip0, initial_a_i, a_max, c1, c2);
%  i=1;
%  for x=-5:0.1:5
%    phia(i)=phi(x);
%    i=i+1;
%  end
%  plot(-5:0.1:5, phia), hold on;
%  plot(a_opt, phi(a_opt), '*');

  
% Initialize with what is given:
a_i_minus_1      = 0;
phi_a_i_minus_1  = phi0;
phip_a_i_minus_1 = phip0;

a_i = initial_a_i;
i   = 1;
while(1)
  phi_a_i = phi_fun(a_i); % Evaluate the function at alpha_i.
  if ((phi_a_i  >  phi0+c1*a_i*phip0) | ...
      ((phi_a_i >= phi_a_i_minus_1) & (i>1)))
   % Already have at-least two distinct points: a=0, a=a_1. 
   %  alpha_opt = zoom(alpha_i_minus_1, alpha_i); 
    a_lo      = a_i_minus_1;
    phi_a_lo  = phi_a_i_minus_1;
    phip_a_lo = phip_a_i_minus_1;
    a_hi      = a_i;
    phi_a_hi  = phi_a_i;
    [a_opt] = zoom (phi_fun, phi_der_fun, phi0, phip0, ... 
                    a_lo, a_hi, phi_a_lo, phi_a_hi, phip_a_lo,  ...
                    c1, c2);           
    return;
  end

  phipai = phi_der_fun(a_i); % phi'(alpha_i)
  if (abs(phipai) <= -c2*phip0)
    a_opt = a_i; 
    return; 
  end
 
  if (phipai >= 0)
    % alpha_opt = zoom(alpha_i, alpha_i_minus_1);     
    a_lo      = a_i;
    phi_a_lo  = phi_a_i;
    phip_a_lo = phipai;
    a_hi      = a_i_minus_1;     % Make sure it is initialized!
    phi_a_hi  = phi_a_i_minus_1;
    [a_opt] = zoom (phi_fun, phi_der_fun, phi0, phip0, ... 
                    a_lo, a_hi, phi_a_lo, phi_a_hi, phip_a_lo,  ...
                    c1, c2);            
    return;            
  end

  % Choose the next trial step based on this one:
  % We must be approachinh a_max so as to terminate the search.
  a_lo      = a_i;
  phi_a_lo  = phi_a_i;
  phip_a_lo = phipai;
  a_hi      = a_i_minus_1;
  phi_a_hi  = phi_a_i_minus_1;
  [a_i_plus_1] = quad_step (a_lo, a_hi, phi_a_lo, phi_a_hi, phip_a_lo);
    
  % If the new a_i is not closer to a_max, force it closer:
  % This will get us closer to a_max.
  if (a_max - a_i_plus_1 > a_max - a_i)
      a_i_plus_1 = (a_max + a_i)/2;
  end      
      
  % Update the previous values:
  a_i_minus_1 = a_i;
  phi_a_i_minus_1  = phi_a_i;
  phip_a_i_minus_1 = phipai;  
  
  % Update the iteration number.
  i = i+1;   
  a_i = a_i_plus_1;
end


