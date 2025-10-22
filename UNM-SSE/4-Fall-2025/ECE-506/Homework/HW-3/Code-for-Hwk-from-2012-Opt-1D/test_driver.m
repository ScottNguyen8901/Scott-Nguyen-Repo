% Clean up Matlab Windows and memory
clear all
close all

% Set the number of trials below:
N = 5;    % Run N random initializations.

% Testing for f1(.): 
f1          = @(x) x^4 - 12*x^3 + 47*x^2 - 60*x;
f1_roots    = [0 3 4 5]; % roots => only test for the first one.
f1_interval = [-5 2];    % assume that we know we are in this interval.

% Testing for f1_der(.):
f1_x        = @(x) 4*x^3 - 36*x^2 + 94*x - 60;
f1_x_zeros  = [0.943 3.456 4.601]; % f1'(x)=0 there.
f1_x_interval = [-5 2]; % assume known.

% Need the second derivative for f1:
f1_xx = @(x) 12*x^2 - 72*x + 94;

% Call for f1, f1_x:
test_zero_funs(f1,   f1_x,  f1_roots,   f1_interval, N);
test_zero_funs(f1_x, f1_xx, f1_x_zeros, f1_x_interval, N);



% Testing for f2(.):
f2          = @(x) x^2 - 3;
f2_roots    = [-sqrt(3) sqrt(3)]; % roots => only test for the first one.
f2_interval = [-4 1]; % Given.

% Testing for f2_x(.):
f2_x      = @(x) 2*x; 
f2_x_zero = [0];        % minimum for x=0.
f2_x_interval = [-1 2]; % Given.

% Need the second derivative for f2(.):
f2_xx = @(x) 2;

% Call for f2, f2_x:
test_zero_funs(f2,   f2_x,  f2_roots,  f2_interval, N);
test_zero_funs(f2_x, f2_xx, f2_x_zero, f2_x_interval, N);



% Testing for f3(.):
f3          = @(x) x^2 - 1;
f3_roots    = [-1 1];
f3_interval = [-5 0];

% Testing for f3_der(.):
f3_x = @(x) 2*x;
f3_x_zero = [0];
f3_x_interval = [-0.9 0.5];

% Need the second derivative for f3:
f3_xx = @(x) 2;

% Call for f3, f3_x:
test_zero_funs(f3,   f3_x,  f3_roots,  f3_interval, N);
test_zero_funs(f3_x, f3_xx, f3_x_zero, f3_x_interval, N);



% Testing for f4(.):
f4          = @(x) x^2 - 2*x + 1;
f4_roots    = [1];
f4_interval = [0 1.5];

% Testing for f4_der(.):
f4_x = @(x) 2*x - 2;
f4_x_zero = [1];
f4_x_interval = [0 1.6];

% Need the second derivative for f4:
f4_xx = @(x) 2;

% Call for f4, f4_x:
test_zero_funs(f4,   f4_x,  f4_roots,  f4_interval, N);
test_zero_funs(f4_x, f4_xx, f4_x_zero, f4_x_interval, N);




% Testing for f5(.):
f5       = @(x) atan(x);
f5_roots = [0];
f5_interval = [-100 500];

% Define the derivative:
f5_x = @(x) 1/(1+x*x);

% Run with a large and a small interval:
test_zero_funs(f5, f5_x, f5_roots, f5_interval, N);
test_zero_funs(f5, f5_x, f5_roots, [-0.8 0.5], N);
