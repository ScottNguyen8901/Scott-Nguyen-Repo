function [min_p] = min_quad (g, B, Delta)
%
% [min_p] = min_quad (g, B, Delta)
% minimizes a quadratic function by evaluating many points.
%
% Input:
%  g:     the gradient
%  B:     the Hessian
%  Delta: Delta for a circular trust region
%
% Output:
%  min_p: contains the point where
%           (g.')*p + 0.5*(p^t)*B*p is minimized.
%
% Example:
% fun   = @(x) x(1).^2 + 5*x(2).^2;
% g_fun = @(x) [2*x(1); 10*x(2)];
% Hessian_fun = @(x) [2, 0; 0 10];
% x = [1; 2];
% g = g_fun(x);
% B = Hessian_fun(x);
% Delta = 5;
% [min_p] = min_quad(g, B, Delta);
% fun(x+min_p)
% g_fun(x+min_p)


% Build the quadratic to minimize:
model_fun = @(p) (g.')*p + 0.5*(p.')*B*p;

N = 100;
min_val = 1.0e100;
for x1=linspace(-Delta, Delta, N)
  for x2=linspace(-Delta, Delta, N)
      % Check for constraints:
      if (x1^2 + x2^2 <= Delta^2)
          p = [x1; x2];
          v = model_fun(p);
          if (v < min_val)
              min_val = v;
              min_p   = [x1; x2];
          end
      end
  end
end  