README FILE
===========
The m-files have been saved under three directories:
line-search, trust-region, and conjugate-gradient.

Each m-file includes a help file. Type:
   >>help filename
to get instructions on how to run each file.



line-search
-----------
The primary file in the line-search directory
is min_rosen.m which calls:
    BFGS(), newton(), and steepest_descent()
to minimize the Rosenbrock function.

The strong Wolfe conditions are implemented in
the line_search(), zoom(), and quad_step() functions.

Results are plotted using plot_trajectory().



trust-region
------------
The primary files are prob_4_1.m and prob_4_3.m
that solve the corresponding problems in the book.
Other files include:
* plot_model() plots local quadratic models
* dogleg()     provides the dogleg method
* min_quad()   provides a brute force method for minimizing
               a quardatic function. It evaluates all points
               on a grid.

* trust_region_CG() implements the trust region method with
               the function CG_Steihaug().

* bisection() is taken from homework #1.


conjugate-gradient
------------------
The primary file is prob_5_1.m that solves the corresponding
problem in the book. The Conjugate Gradient method is implemented
in conj_grad.m.
