close all;
clear; 
clc; 

%% Reading Data

% Define the data sets for Pre-Course and Mid-Course
data_sets = {'PC', 'MC'};
groups = {'GA', 'GB'};
modules = 1:4;

% Define file paths and sheet names
group_A_files = {
    'Group A\GA [1.0] Introduction and Rocket Hardware.xlsx', 
    'Group A\GA [2.0] Fundamentals of Rocketry.xlsx',
    'Group A\GA [3.0] Modeling Rocket Mechanics.xlsx', 
    'Group A\GA [4.0] Analysis.xlsx'};

group_B_files = {
    'Group B\GB [1.0] Introduction and Rocket Hardware.xlsx', 
    'Group B\GB [2.0] Fundamentals of Rocketry.xlsx',
    'Group B\GB [3.0] Modeling Rocket Mechanics.xlsx', 
    'Group B\GB [4.0] Analysis.xlsx'};

% Sheet names
sheet_names = {'PC', 'MC'};

% Load Group A and Group B data for both Pre-Course and Mid-Course
for i = 1:4
    % Pre-Course Data
    eval(['data_PC_GA_mod_' num2str(i) ' = readmatrix(group_A_files{' num2str(i) '}, ''Sheet'', sheet_names{1});']);
    eval(['data_PC_GB_mod_' num2str(i) ' = readmatrix(group_B_files{' num2str(i) '}, ''Sheet'', sheet_names{1});']);
    
    % Mid-Course Data
    eval(['data_MC_GA_mod_' num2str(i) ' = readmatrix(group_A_files{' num2str(i) '}, ''Sheet'', sheet_names{2});']);
    eval(['data_MC_GB_mod_' num2str(i) ' = readmatrix(group_B_files{' num2str(i) '}, ''Sheet'', sheet_names{2});']);
end

% Get the number of students for Group A and Group B
num_stud_GA = length(data_PC_GA_mod_1);
num_stud_GB = length(data_PC_GB_mod_1);

%% Solutions for Modules
sols = {
    [4, 4, 1, 2, 2, 4, 1, 1100]', % Module 1
    [2, 3, 1, 3, 3, 1]',          % Module 2
    [2, 11010, 2, 6, 40, 100]',   % Module 3
    [11011, 2, 4, 1, 1]'          % Module 4
};

% Create solution matrices for Group A and Group B for each module
for i = 1:4
    eval(['sol_mat_GA_mod_' num2str(i) ' = repmat(sols{' num2str(i) '}, 1, num_stud_GA)'';']);
    eval(['sol_mat_GB_mod_' num2str(i) ' = repmat(sols{' num2str(i) '}, 1, num_stud_GB)'';']);
end

%% Differences in Scores

% Calculate differences in scores for each combination
for d = 1:length(data_sets)
    for g = 1:length(groups)
        for m = modules
            eval(['diff_' data_sets{d} '_' groups{g} '_mod_' num2str(m) ' = abs(sol_mat_' groups{g} '_mod_' num2str(m) ' - data_' data_sets{d} '_' groups{g} '_mod_' num2str(m) ');']);
        end
    end
end

%% Calculating Scores

for g = 1:length(groups)
    for m = modules
        eval(sprintf('score_PC_%s_mod_%d = 100 * sum(diff_PC_%s_mod_%d == 0, 2) / size(diff_PC_%s_mod_%d, 2);', ...
            groups{g}, m, groups{g}, m, groups{g}, m));
    end
end

% Mid-Course
for g = 1:length(groups)
    for m = modules
        eval(sprintf('score_MC_%s_mod_%d = 100 * sum(diff_MC_%s_mod_%d == 0, 2) / size(diff_MC_%s_mod_%d, 2);', ...
            groups{g}, m, groups{g}, m, groups{g}, m));
    end
end

%% Statistical Analysis

% Calculate MC to PC differences
for g = 1:length(groups)
    for m = modules
        eval(sprintf('MC_to_PC_%s_mod_%d = score_MC_%s_mod_%d - score_PC_%s_mod_%d;', ...
            groups{g}, m, groups{g}, m, groups{g}, m));
    end
end

% Calculate mean scores for Pre-Course and Mid-Course
for g = 1:length(groups)
    for m = modules
        eval(sprintf('mean_PC_%s_mod_%d = mean(score_PC_%s_mod_%d, ''all'');', ...
            groups{g}, m, groups{g}, m));
        eval(sprintf('mean_MC_%s_mod_%d = mean(score_MC_%s_mod_%d, ''all'');', ...
            groups{g}, m, groups{g}, m));
    end
end

% Combine mean values into matrices
for g = 1:length(groups)
    eval(sprintf('mean_val_%s_mod = [mean_PC_%s_mod_1, mean_MC_%s_mod_1; mean_PC_%s_mod_2, mean_MC_%s_mod_2; mean_PC_%s_mod_3, mean_MC_%s_mod_3; mean_PC_%s_mod_4, mean_MC_%s_mod_4];', ...
        groups{g}, groups{g}, groups{g}, groups{g}, groups{g}, groups{g}, groups{g}, groups{g}, groups{g}));
end

% Calculate mean scores across all modules
for g = 1:length(groups)
    eval(sprintf('mean_score_%s_mod = [mean(score_PC_%s_mod_1, 2), mean(score_MC_%s_mod_1, 2); mean(score_PC_%s_mod_2, 2), mean(score_MC_%s_mod_2, 2); mean(score_PC_%s_mod_3, 2), mean(score_MC_%s_mod_3, 2); mean(score_PC_%s_mod_4, 2), mean(score_MC_%s_mod_4, 2)];', ...
        groups{g}, groups{g}, groups{g}, groups{g}, groups{g}, groups{g}, groups{g}, groups{g}, groups{g}));
end

% Calculate mean differences (MC - PC)
for g = 1:length(groups)
    eval(sprintf('mean_MC_to_PC_%s_mod = [mean_MC_%s_mod_1 - mean_PC_%s_mod_1, mean_MC_%s_mod_2 - mean_PC_%s_mod_2, mean_MC_%s_mod_3 - mean_PC_%s_mod_3, mean_MC_%s_mod_4 - mean_PC_%s_mod_4];', ...
        groups{g}, groups{g}, groups{g}, groups{g}, groups{g}, groups{g}, groups{g}, groups{g}, groups{g}));
end


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

% Combine all scores for Group A and B
all_scores_PC_GA = [score_PC_GA_mod_1, score_PC_GA_mod_2, score_PC_GA_mod_3, score_PC_GA_mod_4];
all_scores_PC_GB = [score_PC_GB_mod_1, score_PC_GB_mod_2, score_PC_GB_mod_3, score_PC_GB_mod_4];
all_scores_MC_GA = [score_MC_GA_mod_1, score_MC_GA_mod_2, score_MC_GA_mod_3, score_MC_GA_mod_4];
all_scores_MC_GB = [score_MC_GB_mod_1, score_MC_GB_mod_2, score_MC_GB_mod_3, score_MC_GB_mod_4];

createBoxplots(all_scores_PC_GA, all_scores_PC_GB, 0.15, 'b', 'r');

createBarCharts(mean_MC_to_PC_GA_mod, mean_MC_to_PC_GB_mod, p_values_GA, p_values_GB, alpha);
