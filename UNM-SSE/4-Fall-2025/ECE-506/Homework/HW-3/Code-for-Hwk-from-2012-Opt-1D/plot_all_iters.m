% Called by test_zero_funs.m to plot the results of all the methods
figure;
plot(bis_evals,     abs(bis_x-root1),     '-',  ... % bisection
     newt_evals,    abs(newt_x-root1),    ':',  ... % Newton
     q_newt_evals,  abs(q_newt_x-root1),  '-.', ... % Quasi-Newton with finite diff delta_x
     q_newt2_evals, abs(q_newt2_x-root1), '--', ... % Quasi-Newton using previous value
     ds_evals,      abs(ds_x-root1),      '+'); 
   
legend('Bisection Method', ...
       'Newton with exact derivatives', ...
       'Quasi-Newton with \delta x finite diff.', ...
       'Quasi-Newton using previous value', ...
       'Dennis and Schnabel Method');
     
xlabel('Number of function Evaluations');
ylabel(sprintf('Distance from the root (r=%5.2f)', root1));
title(sprintf('%s', func2str(fun)));