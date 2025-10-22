function [x] = trust_region_CG (fun, G_fun, H_fun, initial_x, initial_Delta, max_Delta)
%
% [x] = trust_region_CG (fun, G_fun, H_fun, initial_x)
% Applies the trust region approach with Conjugate Gradient, to estimate
% a value for x, where fun is minimized.
%
% Inputs:
%  fun:       Original function
%  G_fun:     Gradient function
%  H_fun:     Hessian function
%  initial_x: initial value for x
%  initial_Delta: initial value for Delta
%  max_Delta:     maximum value allowed for Delta
%
% Outputs:
%  x:         estimate of the point that minimizes fun.
%
% Example:
% vec_eps = 1.0e-4;
% fun   = @(x) x(1).^2 + 5*x(2).^2;
% G_fun = @(x) [2*x(1); 10*x(2)];
% H_fun = @(x) [2, 0; 0 10];
% initial_x = [-2; -3];
% initial_Delta = 1;
% max_Delta     = 10;
% [x] = trust_region_CG (fun, G_fun, H_fun, initial_x, initial_Delta, max_Delta)
% disp('Function at x');
% fun(x)
% disp('Gradient at x');
% g_fun(x)


print_trace = 1;
max_iters   = 100;

eta       = 1/4;
vec_eps   = 0.001; % Epsilon for vector length estimation.
Delta     = initial_Delta;

x_k = initial_x;
for k=0:1:max_iters
  % Build the quadratic model for the CG method:      
  B = H_fun(x_k);
  g = G_fun(x_k);
  if (vec_norm(g)<vec_eps) % Exit if at zero-norm
    x = x_k;
    return
  end
  
  [p_k] = CG_Steihaug (vec_eps, g, B, Delta);
  
  % Evaluate the reduction ratio:
  f_k = fun(x_k);
  f_k_plus_1 = fun(x_k + p_k);
  m_k_0      = f_k;
  m_k_pk     = f_k + (g.')*p_k + 0.5*(p_k.')*B*p_k;  
  reduction_ratio = (f_k - f_k_plus_1)/(m_k_0 - m_k_pk);
  
  % Adjust the trust region size if needed.
  if (reduction_ratio < 1/4)      
      % Very small reduction (or increase!)
      Delta = 0.25*Delta;          
  else    
     % Increase the trust region size.
     if (reduction_ratio>3/4) & (abs(vec_norm(p_k)-Delta)<vec_eps)
        Delta = min(2*Delta, max_Delta);         
     end 
  end
        
  % Decide on whether to take the step or not.
  if (reduction_ratio>eta)
      x_k_plus_1 = x_k + p_k;
  else    
      x_k_plus_1 = x_k;
  end    
  
  if (print_trace==1)
    disp(sprintf('Iteration=%d', k));
    disp(sprintf('Reduction ratio=%f, Delta=%f', reduction_ratio, Delta));
    x_k_plus_1
  end
  
  % Update the step information:
  x_k = x_k_plus_1;  
end    
x = x_k;
