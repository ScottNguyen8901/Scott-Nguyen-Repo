This file describes the files associated with solving
non-linear root-finding and optimization problems.

1. Main files:
The list is as follows:
    test_driver.m      main driver file with 5 examples.
    test_zero_funs.m   runs all methods (see below)
    plot_all_iters.m   plots the results from all methods.

2. Method files:    
    bisection.m        bisection method
    newton.m           Newton with exact derivatives
    quasi_newton.m     Newton with finite-differencing at each x
    quasi_newton2.m    Newton with secant method from previous iterate
    ds_method.m        Dennis & Schnabel's method based on quasi_newton2.

Solution Discussion:
The following test cases were considered in test_driver.m:
    f1   = @(x) x^4 - 12*x^3 + 47*x^2 - 60*x;
    f1_x = @(x) 4*x^3 - 36*x^2 + 94*x - 60;
    f2   = @(x) x^2 - 3;
    f2_x = @(x) 2*x; 
    f3   = @(x) x^2 - 1;
    f3_x = @(x) 2*x;
    f4   = @(x) x^2 - 2*x + 1;
    f4_x = @(x) 2*x - 2;
    f5   = @(x) atan(x);

Convergence is demonstrated in all cases except for:
   (i) arctangent(.): 
       Newton-methods for large deviation from the root, they fail.
       On other hand, we also show that for small deviations from
       the root, convergence is guaranteed.
   (ii) f(x)=x^2-2*x+1 has a double root at 1 and does not alternate
       signs between the roots. Thus, the bisection method fails.

In general, for Newton-based methods we see quadratic convergence,
while the distance is generally halved with each evaluation with
the bisection method. 

