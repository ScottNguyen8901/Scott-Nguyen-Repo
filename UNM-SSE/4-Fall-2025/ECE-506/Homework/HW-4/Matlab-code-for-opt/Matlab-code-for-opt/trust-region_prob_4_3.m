clear all
close all

initial_Delta = 1;
max_Delta     = 10;

fun   = @fun_4_3;
G_fun = @Gradient_fun; 
H_fun = @Hessian_fun;

N = 50;
for k=1:10
  initial_x = randn(N,1);
  [x] = trust_region_CG (fun, G_fun, H_fun, initial_x, initial_Delta, max_Delta);

  disp('Function at x');
  fun(x)
  disp('Gradient at x');
  G_fun(x)
end