clear; clc; close all;

%% Reading Data

% Group A Pre-Course Data
data_PC_GA_mod_1 = readmatrix('Group A\GA [1.0] Introduction and Rocket Hardware.xlsx', 'Sheet', 'PC');
data_PC_GA_mod_2 = readmatrix('Group A\GA [2.0] Fundamentals of Rocketry.xlsx', 'Sheet', 'PC');
data_PC_GA_mod_3 = readmatrix('Group A\GA [3.0] Modeling Rocket Mechanics.xlsx', 'Sheet', 'PC');
data_PC_GA_mod_4 = readmatrix('Group A\GA [4.0] Analysis.xlsx', 'Sheet', 'PC');

% Group B Pre-Course Data
data_PC_GB_mod_1 = readmatrix('Group B\GB [1.0] Introduction and Rocket Hardware.xlsx', 'Sheet', 'PC');
data_PC_GB_mod_2 = readmatrix('Group B\GB [2.0] Fundamentals of Rocketry.xlsx', 'Sheet', 'PC');
data_PC_GB_mod_3 = readmatrix('Group B\GB [3.0] Modeling Rocket Mechanics.xlsx', 'Sheet', 'PC');
data_PC_GB_mod_4 = readmatrix('Group B\GB [4.0] Analysis.xlsx', 'Sheet', 'PC');

% Group A Mid-Course Data
data_MC_GA_mod_1 = readmatrix('Group A\GA [1.0] Introduction and Rocket Hardware.xlsx', 'Sheet', 'MC');
data_MC_GA_mod_2 = readmatrix('Group A\GA [2.0] Fundamentals of Rocketry.xlsx', 'Sheet', 'MC');
data_MC_GA_mod_3 = readmatrix('Group A\GA [3.0] Modeling Rocket Mechanics.xlsx', 'Sheet', 'MC');
data_MC_GA_mod_4 = readmatrix('Group A\GA [4.0] Analysis.xlsx', 'Sheet', 'MC');

% Group B Mid-Course Data
data_MC_GB_mod_1 = readmatrix('Group B\GB [1.0] Introduction and Rocket Hardware.xlsx', 'Sheet', 'MC');
data_MC_GB_mod_2 = readmatrix('Group B\GB [2.0] Fundamentals of Rocketry.xlsx', 'Sheet', 'MC');
data_MC_GB_mod_3 = readmatrix('Group B\GB [3.0] Modeling Rocket Mechanics.xlsx', 'Sheet', 'MC');
data_MC_GB_mod_4 = readmatrix('Group B\GB [4.0] Analysis.xlsx', 'Sheet', 'MC');

num_stud_GA = length(data_PC_GA_mod_1);
num_stud_GB = length(data_PC_GB_mod_1);

%% Solutions for Modules
sol_mod_1 = [4, 4, 1, 2, 2, 4, 1, 1100]';
sol_mod_2 = [2, 3, 1, 3, 3, 1]';
sol_mod_3 = [2, 11010, 2, 6, 40, 100]';
sol_mod_4 = [11011, 2, 4, 1, 1]';

% Group A solution matrix for modules 1-4
sol_mat_GA_mod_1 = repmat(sol_mod_1, 1, num_stud_GA)';
sol_mat_GA_mod_2 = repmat(sol_mod_2, 1, num_stud_GA)';
sol_mat_GA_mod_3 = repmat(sol_mod_3, 1, num_stud_GA)';
sol_mat_GA_mod_4 = repmat(sol_mod_4, 1, num_stud_GA)';

% Group B solution matrix for modules 1-4
sol_mat_GB_mod_1 = repmat(sol_mod_1, 1, num_stud_GB)';
sol_mat_GB_mod_2 = repmat(sol_mod_2, 1, num_stud_GB)';
sol_mat_GB_mod_3 = repmat(sol_mod_3, 1, num_stud_GB)';
sol_mat_GB_mod_4 = repmat(sol_mod_4, 1, num_stud_GB)';

%% Differences in Scores
% Pre-Course
diff_PC_GA_mod_1 = abs(sol_mat_GA_mod_1 - data_PC_GA_mod_1);
diff_PC_GA_mod_2 = abs(sol_mat_GA_mod_2 - data_PC_GA_mod_2);
diff_PC_GA_mod_3 = abs(sol_mat_GA_mod_3 - data_PC_GA_mod_3);
diff_PC_GA_mod_4 = abs(sol_mat_GA_mod_4 - data_PC_GA_mod_4);

diff_PC_GB_mod_1 = abs(sol_mat_GB_mod_1 - data_PC_GB_mod_1);
diff_PC_GB_mod_2 = abs(sol_mat_GB_mod_2 - data_PC_GB_mod_2);
diff_PC_GB_mod_3 = abs(sol_mat_GB_mod_3 - data_PC_GB_mod_3);
diff_PC_GB_mod_4 = abs(sol_mat_GB_mod_4 - data_PC_GB_mod_4);

% Mid-Course
diff_MC_GA_mod_1 = abs(sol_mat_GA_mod_1 - data_MC_GA_mod_1);
diff_MC_GA_mod_2 = abs(sol_mat_GA_mod_2 - data_MC_GA_mod_2);
diff_MC_GA_mod_3 = abs(sol_mat_GA_mod_3 - data_MC_GA_mod_3);
diff_MC_GA_mod_4 = abs(sol_mat_GA_mod_4 - data_MC_GA_mod_4);

diff_MC_GB_mod_1 = abs(sol_mat_GB_mod_1 - data_MC_GB_mod_1);
diff_MC_GB_mod_2 = abs(sol_mat_GB_mod_2 - data_MC_GB_mod_2);
diff_MC_GB_mod_3 = abs(sol_mat_GB_mod_3 - data_MC_GB_mod_3);
diff_MC_GB_mod_4 = abs(sol_mat_GB_mod_4 - data_MC_GB_mod_4);

%% Calculating Scores
% Pre-Course
score_PC_GA_mod_1 = 100*sum(diff_PC_GA_mod_1 == 0, 2)/size(diff_PC_GA_mod_1, 2);
score_PC_GA_mod_2 = 100*sum(diff_PC_GA_mod_2 == 0, 2)/size(diff_PC_GA_mod_2, 2);
score_PC_GA_mod_3 = 100*sum(diff_PC_GA_mod_3 == 0, 2)/size(diff_PC_GA_mod_3, 2);
score_PC_GA_mod_4 = 100*sum(diff_PC_GA_mod_4 == 0, 2)/size(diff_PC_GA_mod_4, 2);

score_PC_GB_mod_1 = 100*sum(diff_PC_GB_mod_1 == 0, 2)/size(diff_PC_GB_mod_1, 2);
score_PC_GB_mod_2 = 100*sum(diff_PC_GB_mod_2 == 0, 2)/size(diff_PC_GB_mod_2, 2);
score_PC_GB_mod_3 = 100*sum(diff_PC_GB_mod_3 == 0, 2)/size(diff_PC_GB_mod_3, 2);
score_PC_GB_mod_4 = 100*sum(diff_PC_GB_mod_4 == 0, 2)/size(diff_PC_GB_mod_4, 2);

% Mid-Course
score_MC_GA_mod_1 = 100*sum(diff_MC_GA_mod_1 == 0, 2)/size(diff_MC_GA_mod_1, 2);
score_MC_GA_mod_2 = 100*sum(diff_MC_GA_mod_2 == 0, 2)/size(diff_MC_GA_mod_2, 2);
score_MC_GA_mod_3 = 100*sum(diff_MC_GA_mod_3 == 0, 2)/size(diff_MC_GA_mod_3, 2);
score_MC_GA_mod_4 = 100*sum(diff_MC_GA_mod_4 == 0, 2)/size(diff_MC_GA_mod_4, 2);

score_MC_GB_mod_1 = 100*sum(diff_MC_GB_mod_1 == 0, 2)/size(diff_MC_GB_mod_1, 2);
score_MC_GB_mod_2 = 100*sum(diff_MC_GB_mod_2 == 0, 2)/size(diff_MC_GB_mod_2, 2);
score_MC_GB_mod_3 = 100*sum(diff_MC_GB_mod_3 == 0, 2)/size(diff_MC_GB_mod_3, 2);
score_MC_GB_mod_4 = 100*sum(diff_MC_GB_mod_4 == 0, 2)/size(diff_MC_GB_mod_4, 2);

%% Statistical Analysis

MC_to_PC_GA_mod_1 = score_MC_GA_mod_1 - score_PC_GA_mod_1;
MC_to_PC_GA_mod_2 = score_MC_GA_mod_2 - score_PC_GA_mod_2;
MC_to_PC_GA_mod_3 = score_MC_GA_mod_3 - score_PC_GA_mod_3;
MC_to_PC_GA_mod_4 = score_MC_GA_mod_4 - score_PC_GA_mod_4;

MC_to_PC_GB_mod_1 = score_MC_GB_mod_1 - score_PC_GB_mod_1;
MC_to_PC_GB_mod_2 = score_MC_GB_mod_2 - score_PC_GB_mod_2;
MC_to_PC_GB_mod_3 = score_MC_GB_mod_3 - score_PC_GB_mod_3;
MC_to_PC_GB_mod_4 = score_MC_GB_mod_4 - score_PC_GB_mod_4;

mean_PC_GA_mod_1 = mean(score_PC_GA_mod_1, 'all');
mean_PC_GA_mod_2 = mean(score_PC_GA_mod_2, 'all');
mean_PC_GA_mod_3 = mean(score_PC_GA_mod_3, 'all');
mean_PC_GA_mod_4 = mean(score_PC_GA_mod_4, 'all');

mean_PC_GB_mod_1 = mean(score_PC_GB_mod_1, 'all');
mean_PC_GB_mod_2 = mean(score_PC_GB_mod_2, 'all');
mean_PC_GB_mod_3 = mean(score_PC_GB_mod_3, 'all');
mean_PC_GB_mod_4 = mean(score_PC_GB_mod_4, 'all');

mean_MC_GA_mod_1 = mean(score_MC_GA_mod_1, 'all');
mean_MC_GA_mod_2 = mean(score_MC_GA_mod_2, 'all');
mean_MC_GA_mod_3 = mean(score_MC_GA_mod_3, 'all');
mean_MC_GA_mod_4 = mean(score_MC_GA_mod_4, 'all');

mean_MC_GB_mod_1 = mean(score_MC_GB_mod_1, 'all');
mean_MC_GB_mod_2 = mean(score_MC_GB_mod_2, 'all');
mean_MC_GB_mod_3 = mean(score_MC_GB_mod_3, 'all');
mean_MC_GB_mod_4 = mean(score_MC_GB_mod_4, 'all');

mean_val_GA_mod = [mean_PC_GA_mod_1, mean_MC_GA_mod_1;
                   mean_PC_GA_mod_2, mean_MC_GA_mod_2;
                   mean_PC_GA_mod_3, mean_MC_GA_mod_3;
                   mean_PC_GA_mod_4, mean_MC_GA_mod_4];
                   
mean_val_GB_mod = [mean_PC_GB_mod_1, mean_MC_GB_mod_1;
                   mean_PC_GB_mod_2, mean_MC_GB_mod_2;
                   mean_PC_GB_mod_3, mean_MC_GB_mod_3;
                   mean_PC_GB_mod_4, mean_MC_GB_mod_4];

mean_score_GA_mod = [mean(score_PC_GA_mod_1, 2), mean(score_MC_GA_mod_1, 2);
                     mean(score_PC_GA_mod_2, 2), mean(score_MC_GA_mod_2, 2);
                     mean(score_PC_GA_mod_3, 2), mean(score_MC_GA_mod_3, 2);
                     mean(score_PC_GA_mod_4, 2), mean(score_MC_GA_mod_4, 2)];

mean_score_GB_mod = [mean(score_PC_GB_mod_1, 2), mean(score_MC_GB_mod_1, 2);
                     mean(score_PC_GB_mod_2, 2), mean(score_MC_GB_mod_2, 2);
                     mean(score_PC_GB_mod_3, 2), mean(score_MC_GB_mod_3, 2);
                     mean(score_PC_GB_mod_4, 2), mean(score_MC_GB_mod_4, 2)];

mean_MC_to_PC_GA_mod = [mean_MC_GA_mod_1 - mean_PC_GA_mod_1, 
                        mean_MC_GA_mod_2 - mean_PC_GA_mod_2, 
                        mean_MC_GA_mod_3 - mean_PC_GA_mod_3, 
                        mean_MC_GA_mod_4 - mean_PC_GA_mod_4];

mean_MC_to_PC_GB_mod = [mean_MC_GB_mod_1 - mean_PC_GB_mod_1, 
                        mean_MC_GB_mod_2 - mean_PC_GB_mod_2, 
                        mean_MC_GB_mod_3 - mean_PC_GB_mod_3, 
                        mean_MC_GB_mod_4 - mean_PC_GB_mod_4];

alpha = 0.01;

% Pre-course comparisons
[H_PC_mod_1, p_PC_mod_1, ~] = ttest2(mean(score_PC_GA_mod_1, 2), mean(score_PC_GB_mod_1, 2), 'Tail', 'both', 'Alpha', alpha);
[H_PC_mod_2, p_PC_mod_2, ~] = ttest2(mean(score_PC_GA_mod_2, 2), mean(score_PC_GB_mod_2, 2), 'Tail', 'both', 'Alpha', alpha);
[H_PC_mod_3, p_PC_mod_3, ~] = ttest2(mean(score_PC_GA_mod_3, 2), mean(score_PC_GB_mod_3, 2), 'Tail', 'both', 'Alpha', alpha);
[H_PC_mod_4, p_PC_mod_4, ~] = ttest2(mean(score_PC_GA_mod_4, 2), mean(score_PC_GB_mod_4, 2), 'Tail', 'both', 'Alpha', alpha);

% Mid-course comparisons
[H_MC_mod_1, p_MC_mod_1, ~] = ttest2(mean(score_MC_GA_mod_1, 2), mean(score_MC_GB_mod_1, 2), 'Tail', 'both', 'Alpha', alpha);
[H_MC_mod_2, p_MC_mod_2, ~] = ttest2(mean(score_MC_GA_mod_2, 2), mean(score_MC_GB_mod_2, 2), 'Tail', 'both', 'Alpha', alpha);
[H_MC_mod_3, p_MC_mod_3, ~] = ttest2(mean(score_MC_GA_mod_3, 2), mean(score_MC_GB_mod_3, 2), 'Tail', 'both', 'Alpha', alpha);
[H_MC_mod_4, p_MC_mod_4, ~] = ttest2(mean(score_MC_GA_mod_4, 2), mean(score_MC_GB_mod_4, 2), 'Tail', 'both', 'Alpha', alpha);

% Pre-course vs Mid-course comparisons for each module and group
% Group A
[H_PC_MC_GA_mod_1, p_PC_MC_GA_mod_1, ~] = ttest2(score_PC_GA_mod_1, score_MC_GA_mod_1, 'Tail', 'both', 'Alpha', alpha);
[H_PC_MC_GA_mod_2, p_PC_MC_GA_mod_2, ~] = ttest2(score_PC_GA_mod_2, score_MC_GA_mod_2, 'Tail', 'both', 'Alpha', alpha);
[H_PC_MC_GA_mod_3, p_PC_MC_GA_mod_3, ~] = ttest2(score_PC_GA_mod_3, score_MC_GA_mod_3, 'Tail', 'both', 'Alpha', alpha);
[H_PC_MC_GA_mod_4, p_PC_MC_GA_mod_4, ~] = ttest2(score_PC_GA_mod_4, score_MC_GA_mod_4, 'Tail', 'both', 'Alpha', alpha);

% Group B
[H_PC_MC_GB_mod_1, p_PC_MC_GB_mod_1, ~] = ttest2(score_PC_GB_mod_1, score_MC_GB_mod_1, 'Tail', 'both', 'Alpha', alpha);
[H_PC_MC_GB_mod_2, p_PC_MC_GB_mod_2, ~] = ttest2(score_PC_GB_mod_2, score_MC_GB_mod_2, 'Tail', 'both', 'Alpha', alpha);
[H_PC_MC_GB_mod_3, p_PC_MC_GB_mod_3, ~] = ttest2(score_PC_GB_mod_3, score_MC_GB_mod_3, 'Tail', 'both', 'Alpha', alpha);
[H_PC_MC_GB_mod_4, p_PC_MC_GB_mod_4, ~] = ttest2(score_PC_GB_mod_4, score_MC_GB_mod_4, 'Tail', 'both', 'Alpha', alpha);

% Store results in vectors
p_values_GA = [p_PC_MC_GA_mod_1, p_PC_MC_GA_mod_2, p_PC_MC_GA_mod_3, p_PC_MC_GA_mod_4];
p_values_GB = [p_PC_MC_GB_mod_1, p_PC_MC_GB_mod_2, p_PC_MC_GB_mod_3, p_PC_MC_GB_mod_4];

H_values_GA = [H_PC_MC_GA_mod_1, H_PC_MC_GA_mod_2, H_PC_MC_GA_mod_3, H_PC_MC_GA_mod_4];
H_values_GB = [H_PC_MC_GB_mod_1, H_PC_MC_GB_mod_2, H_PC_MC_GB_mod_3, H_PC_MC_GB_mod_4];

%% Printing 
% Store the results in a MATLAB table
results_table = {'Module', 'P-Value', 'Reject Null'};
results_table{2,1} = 'Module 1';
results_table{3,1} = 'Module 2';
results_table{4,1} = 'Module 3';
results_table{5,1} = 'Module 4';

% Fill in the table with data
results_table{2,2} = p_PC_mod_1;
results_table{3,2} = p_PC_mod_2;
results_table{4,2} = p_PC_mod_3;
results_table{5,2} = p_PC_mod_4;

results_table{2,3} = ternary(H_PC_mod_1, 'Yes', 'No');
results_table{3,3} = ternary(H_PC_mod_2, 'Yes', 'No');
results_table{4,3} = ternary(H_PC_mod_3, 'Yes', 'No');
results_table{5,3} = ternary(H_PC_mod_4, 'Yes', 'No');

% Print the table
fprintf('Pre-Course Module Comparison Results:\n');
fprintf('%-10s %-15s %-15s\n', results_table{1,:});
for i = 2:size(results_table,1)
    fprintf('%-10s %-15.4f %-15s\n', results_table{i,1}, results_table{i,2}, results_table{i,3});
end

% Store the results in a MATLAB table
results_table_mid = {'Module', 'P-Value', 'Reject Null'};
results_table_mid{2,1} = 'Module 1';
results_table_mid{3,1} = 'Module 2';
results_table_mid{4,1} = 'Module 3';
results_table_mid{5,1} = 'Module 4';

% Fill in the table with data
results_table_mid{2,2} = p_MC_mod_1;
results_table_mid{3,2} = p_MC_mod_2;
results_table_mid{4,2} = p_MC_mod_3;
results_table_mid{5,2} = p_MC_mod_4;

results_table_mid{2,3} = ternary(H_MC_mod_1, 'Yes', 'No');
results_table_mid{3,3} = ternary(H_MC_mod_2, 'Yes', 'No');
results_table_mid{4,3} = ternary(H_MC_mod_3, 'Yes', 'No');
results_table_mid{5,3} = ternary(H_MC_mod_4, 'Yes', 'No');

% Print the table
fprintf('\nMid-Course Module Comparison Results:\n');
fprintf('%-10s %-15s %-15s\n', results_table_mid{1,:});
for i = 2:size(results_table_mid,1)
    fprintf('%-10s %-15.4f %-15s\n', results_table_mid{i,1}, results_table_mid{i,2}, results_table_mid{i,3});
end

% Store the results in a MATLAB table
results_table_GA = {'Module', 'Mean PC', 'Mean MC', 'Mean Change'};
results_table_GA{2,1} = 'Module 1';
results_table_GA{3,1} = 'Module 2';
results_table_GA{4,1} = 'Module 3';
results_table_GA{5,1} = 'Module 4';

% Fill in the table with data
results_table_GA{2,2} = mean(score_PC_GA_mod_1, 'all');
results_table_GA{3,2} = mean(score_PC_GA_mod_2, 'all');
results_table_GA{4,2} = mean(score_PC_GA_mod_3, 'all');
results_table_GA{5,2} = mean(score_PC_GA_mod_4, 'all');

results_table_GA{2,3} = mean(score_MC_GA_mod_1, 'all');
results_table_GA{3,3} = mean(score_MC_GA_mod_2, 'all');
results_table_GA{4,3} = mean(score_MC_GA_mod_3, 'all');
results_table_GA{5,3} = mean(score_MC_GA_mod_4, 'all');

results_table_GA{2,4} = mean_MC_to_PC_GA_mod(1);
results_table_GA{3,4} = mean_MC_to_PC_GA_mod(2);
results_table_GA{4,4} = mean_MC_to_PC_GA_mod(3);
results_table_GA{5,4} = mean_MC_to_PC_GA_mod(4);

% Print the table
fprintf('\nGroup A Mean Scores and Changes:\n');
fprintf('%-10s %-10s %-10s %-10s\n', results_table_GA{1,:});
for i = 2:size(results_table_GA,1)
    fprintf('%-10s %-10.4f %-10.4f %-10.4f\n', results_table_GA{i,:});
end

% Store the results in a MATLAB table
results_table_GB = {'Module', 'Mean PC', 'Mean MC', 'Mean Change'};
results_table_GB{2,1} = 'Module 1';
results_table_GB{3,1} = 'Module 2';
results_table_GB{4,1} = 'Module 3';
results_table_GB{5,1} = 'Module 4';

% Fill in the table with data
results_table_GB{2,2} = mean(score_PC_GB_mod_1, 'all');
results_table_GB{3,2} = mean(score_PC_GB_mod_2, 'all');
results_table_GB{4,2} = mean(score_PC_GB_mod_3, 'all');
results_table_GB{5,2} = mean(score_PC_GB_mod_4, 'all');

results_table_GB{2,3} = mean(score_MC_GB_mod_1, 'all');
results_table_GB{3,3} = mean(score_MC_GB_mod_2, 'all');
results_table_GB{4,3} = mean(score_MC_GB_mod_3, 'all');
results_table_GB{5,3} = mean(score_MC_GB_mod_4, 'all');

results_table_GB{2,4} = mean_MC_to_PC_GB_mod(1);
results_table_GB{3,4} = mean_MC_to_PC_GB_mod(2);
results_table_GB{4,4} = mean_MC_to_PC_GB_mod(3);
results_table_GB{5,4} = mean_MC_to_PC_GB_mod(4);

% Print the table
fprintf('\nGroup B Mean Scores and Changes:\n');
fprintf('%-10s %-10s %-10s %-10s\n', results_table_GB{1,:});
for i = 2:size(results_table_GB,1)
    fprintf('%-10s %-10.4f %-10.4f %-10.4f\n', results_table_GB{i,:});
end

% Store the results in a MATLAB table for Group A
results_table_GA = {'Module', 'P-Value', 'Reject Null'};
results_table_GA{2,1} = 'Module 1';
results_table_GA{3,1} = 'Module 2';
results_table_GA{4,1} = 'Module 3';
results_table_GA{5,1} = 'Module 4';

% Fill in the table with data for Group A
results_table_GA{2,2} = p_PC_MC_GA_mod_1;
results_table_GA{3,2} = p_PC_MC_GA_mod_2;
results_table_GA{4,2} = p_PC_MC_GA_mod_3;
results_table_GA{5,2} = p_PC_MC_GA_mod_4;

results_table_GA{2,3} = ternary(H_PC_MC_GA_mod_1, 'Yes', 'No');
results_table_GA{3,3} = ternary(H_PC_MC_GA_mod_2, 'Yes', 'No');
results_table_GA{4,3} = ternary(H_PC_MC_GA_mod_3, 'Yes', 'No');
results_table_GA{5,3} = ternary(H_PC_MC_GA_mod_4, 'Yes', 'No');

% Print the table for Group A
fprintf('Group A Pre-Course vs Mid-Course Module Comparison Results:\n');
fprintf('%-10s %-15s %-15s\n', results_table_GA{1,:});
for i = 2:size(results_table_GA,1)
    fprintf('%-10s %-15.1e %-15s\n', results_table_GA{i,1}, results_table_GA{i,2}, results_table_GA{i,3});
end

% Store the results in a MATLAB table for Group B
results_table_GB = {'Module', 'P-Value', 'Reject Null'};
results_table_GB{2,1} = 'Module 1';
results_table_GB{3,1} = 'Module 2';
results_table_GB{4,1} = 'Module 3';
results_table_GB{5,1} = 'Module 4';

% Fill in the table with data for Group B
results_table_GB{2,2} = p_PC_MC_GB_mod_1;
results_table_GB{3,2} = p_PC_MC_GB_mod_2;
results_table_GB{4,2} = p_PC_MC_GB_mod_3;
results_table_GB{5,2} = p_PC_MC_GB_mod_4;

results_table_GB{2,3} = ternary(H_PC_MC_GB_mod_1, 'Yes', 'No');
results_table_GB{3,3} = ternary(H_PC_MC_GB_mod_2, 'Yes', 'No');
results_table_GB{4,3} = ternary(H_PC_MC_GB_mod_3, 'Yes', 'No');
results_table_GB{5,3} = ternary(H_PC_MC_GB_mod_4, 'Yes', 'No');

% Print the table for Group B
fprintf('Group B Pre-Course vs Mid-Course Module Comparison Results:\n');
fprintf('%-10s %-15s %-15s\n', results_table_GB{1,:});
for i = 2:size(results_table_GB,1)
    fprintf('%-10s %-15.1e %-15s\n', results_table_GB{i,1}, results_table_GB{i,2}, results_table_GB{i,3});
end

%% Plotting

fs = 16; % Define the font size

edges = 0:10:100; % Define edges for histogram bins

for module = 1:4
    % Pre-course scores subplot
    subplot(4, 2, 2*(module-1) + 1);
    hold on;
    histogram(eval(['score_PC_GA_mod_', num2str(module)]), edges, 'FaceAlpha', 0.5, 'FaceColor', 'b');
    histogram(eval(['score_PC_GB_mod_', num2str(module)]), edges, 'FaceAlpha', 0.5, 'FaceColor', 'r');
    if module == 4 % Only for the last row
        xlabel('Scores', 'FontName', 'Times New Roman', 'FontSize', fs);
    end
    ylabel(['Module ', num2str(module)], 'FontName', 'Times New Roman', 'FontSize', fs);

    grid on;
    set(gca, 'FontName', 'Times New Roman', 'FontSize', fs); % Set font for axes

    % Mid-course scores subplot
    subplot(4, 2, 2*module);
    hold on;
    histogram(eval(['score_MC_GA_mod_', num2str(module)]), edges, 'FaceAlpha', 0.5, 'FaceColor', 'b');
    histogram(eval(['score_MC_GB_mod_', num2str(module)]), edges, 'FaceAlpha', 0.5, 'FaceColor', 'r');
    if module == 4 % Only for the last row
        xlabel('Scores', 'FontName', 'Times New Roman', 'FontSize', fs);
    end

    grid on;
    set(gca, 'FontName', 'Times New Roman', 'FontSize', fs); % Set font for axes
end

% Combine all scores into a matrix for Group A and B
all_scores_PC_GA = [score_PC_GA_mod_1, score_PC_GA_mod_2, score_PC_GA_mod_3, score_PC_GA_mod_4];
all_scores_PC_GB = [score_PC_GB_mod_1, score_PC_GB_mod_2, score_PC_GB_mod_3, score_PC_GB_mod_4];

% Create boxplots for all four modules for Group A with blue color
figure;
subplot(1,2,1)
pos_1 = [1, 2, 3, 4] - 0.15; % Positions for Group A
h_GA = boxplot(all_scores_PC_GA, 'Labels', {'Module 1', 'Module 2', 'Module 3', 'Module 4'}, 'Colors', 'b', 'Widths', 0.3, 'Positions', pos_1);
hold on;

% Create boxplots for all four modules for Group B with red color
pos_2 = pos_1 + 0.3; % Positions for Group B
h_GB = boxplot(all_scores_PC_GB, 'Labels', {'Module 1', 'Module 2', 'Module 3', 'Module 4'}, 'Colors', 'r', 'Widths', 0.3, 'Positions', pos_2);

ylabel('Scores (%)', 'FontName', 'Times New Roman', 'FontSize', 24);

% Example of adjusting XTickLabel location
x_label = {'Module 1', 'Module 2', 'Module 3', 'Module 4'};
x_tick_values = 1:numel(x_label);
set(gca, 'XTickLabel', x_label, 'XTick', x_tick_values, 'FontSize', 24, 'FontName', 'Times New Roman');

% Modifying line widths
set(findobj(h_GA, 'type', 'line'), 'LineWidth', 2.5); % Adjust the line width as needed
set(findobj(h_GB, 'type', 'line'), 'LineWidth', 2.5); % Adjust the line width as needed

% Set the color of outliers to blue for Group A
outliers_GA = findobj(h_GA, 'Tag', 'Outliers');
set(outliers_GA, 'MarkerEdgeColor', 'b');

% Set the color of outliers to red for Group B
outliers_GB = findobj(h_GB, 'Tag', 'Outliers');
set(outliers_GB, 'MarkerEdgeColor', 'r');

% Add mean values to the plot for Group A
for i = 1:numel(pos_1)
    text(pos_1(i), mean(all_scores_PC_GA(:,i)) + 5, sprintf('%.2f', mean(all_scores_PC_GA(:,i))), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10, 'Color', 'b', 'FontName', 'Times New Roman');
end

% Add mean values to the plot for Group B
for i = 1:numel(pos_2)
    text(pos_2(i), mean(all_scores_PC_GB(:,i)) + 5, sprintf('%.2f', mean(all_scores_PC_GB(:,i))), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10, 'Color', 'r', 'FontName', 'Times New Roman');
end

axis square
grid on;

% Combine all scores into a matrix for Group A and B
all_scores_MC_GA = [score_MC_GA_mod_1, score_MC_GA_mod_2, score_MC_GA_mod_3, score_MC_GA_mod_4];
all_scores_MC_GB = [score_MC_GB_mod_1, score_MC_GB_mod_2, score_MC_GB_mod_3, score_MC_GB_mod_4];

% Create boxplots for all four modules for Group A with blue color
subplot(1,2,2)
pos_1 = [1, 2, 3, 4] - 0.15; % Positions for Group A
h_GA = boxplot(all_scores_MC_GA, 'Labels', {'Module 1', 'Module 2', 'Module 3', 'Module 4'}, 'Colors', 'b', 'Widths', 0.3, 'Positions', pos_1);
hold on;

% Create boxplots for all four modules for Group B with red color
pos_2 = pos_1 + 0.3; % Positions for Group B
h_GB = boxplot(all_scores_MC_GB, 'Labels', {'Module 1', 'Module 2', 'Module 3', 'Module 4'}, 'Colors', 'r', 'Widths', 0.3, 'Positions', pos_2);

ylabel('Scores (%)', 'FontName', 'Times New Roman', 'FontSize', 24);

% Example of adjusting XTickLabel location
x_label = {'Module 1', 'Module 2', 'Module 3', 'Module 4'};
x_tick_values = 1:numel(x_label);
set(gca, 'XTickLabel', x_label, 'XTick', x_tick_values, 'FontSize', 24, 'FontName', 'Times New Roman');

% Modifying line widths
set(findobj(h_GA, 'type', 'line'), 'LineWidth', 2.5); % Adjust the line width as needed
set(findobj(h_GB, 'type', 'line'), 'LineWidth', 2.5); % Adjust the line width as needed

% Set the color of outliers to blue for Group A
outliers_GA = findobj(h_GA, 'Tag', 'Outliers');
set(outliers_GA, 'MarkerEdgeColor', 'b');

% Set the color of outliers to red for Group B
outliers_GB = findobj(h_GB, 'Tag', 'Outliers');
set(outliers_GB, 'MarkerEdgeColor', 'r');

% Add mean values to the plot for Group A
for i = 1:numel(pos_1)
    text(pos_1(i), mean(all_scores_MC_GA(:,i)) + 5, sprintf('%.2f', mean(all_scores_MC_GA(:,i))), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10, 'Color', 'b', 'FontName', 'Times New Roman');
end

% Add mean values to the plot for Group B
for i = 1:numel(pos_2)
    text(pos_2(i), mean(all_scores_MC_GB(:,i)) + 3, sprintf('%.2f', mean(all_scores_MC_GB(:,i))), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10, 'Color', 'r', 'FontName', 'Times New Roman');
end

% Add legend
legend([h_GA(1), h_GB(1)], {'Group A', 'Group B'}, 'Location', 'best', 'FontSize', 16, 'FontName', 'Times New Roman');

axis square
grid on;

% Plotting bar charts overlaying Group A and Group B scores for each module
modules = 1:4; % Number of modules

figure;
hold on;

% Set font properties
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);

bar(modules - 0.2, mean_MC_to_PC_GA_mod, 0.4, 'b'); % Group A bars
bar(modules + 0.2, mean_MC_to_PC_GB_mod, 0.4, 'r'); % Group B bars

% Add bar values on top of the bars with the Greek letter delta
for i = 1:length(modules)
    text(modules(i) -0.2, mean_MC_to_PC_GA_mod(i) + 0.1, sprintf('%.2f', mean_MC_to_PC_GA_mod(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 16, 'Color', 'b', 'FontName', 'Times New Roman', 'Interpreter', 'tex');
    text(modules(i) +0.2, mean_MC_to_PC_GB_mod(i) + 0.1, sprintf('%.2f', mean_MC_to_PC_GB_mod(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 16, 'Color', 'r', 'FontName', 'Times New Roman', 'Interpreter', 'tex');

    text(modules(i)-0.25, -2, ...
        sprintf('%.1e\n%s', p_values_GA(i), decision(p_values_GA(i), alpha)), ...
        'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', 'b', 'FontName', 'Times New Roman');

    text(modules(i)+0.25,- 2, ...
        sprintf('%.1e\n%s', p_values_GB(i), decision(p_values_GB(i), alpha)), ...
        'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', 'r', 'FontName', 'Times New Roman');
end

text(0.25, -1.5, ...
    sprintf('P-Val:'), 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', 'k', 'FontName', 'Times New Roman');

text(0.25, -3, ...
    sprintf('Rej H_0:'), 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', 'k', 'FontName', 'Times New Roman');

xlim([-0.25 5])
ylim([-5 33])
ylabel('Change in Score (%)', 'FontName', 'Times New Roman', 'FontSize', 24);
legend('Group A', 'Group B', 'FontName', 'Times New Roman', 'FontSize', 24);
xticks(modules);
xticklabels({'Module 1', 'Module 2', 'Module 3', 'Module 4'});
set(gca, 'XTickLabel', get(gca, 'XTickLabel'), 'FontName', 'Times New Roman', 'FontSize', 24);
grid on;
axis square
hold off;

