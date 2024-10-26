clear; clc; close all;

% Group A Module 1-4 Data
data_PC_GA_mod_1 = readmatrix('G:\My Drive\UIUC\Spring 2024\DodNDEP\Results\Data Analysis\Group A\GA [1.0] Introduction and Rocket Hardware.xlsx', 'Sheet', 'MC');
data_PC_GA_mod_2 = readmatrix('G:\My Drive\UIUC\Spring 2024\DodNDEP\Results\Data Analysis\Group A\GA [2.0] Fundamentals of Rocketry.xlsx', 'Sheet', 'MC');
data_PC_GA_mod_3 = readmatrix('G:\My Drive\UIUC\Spring 2024\DodNDEP\Results\Data Analysis\Group A\GA [3.0] Modeling Rocket Mechanics.xlsx', 'Sheet', 'MC');
data_PC_GA_mod_4 = readmatrix('G:\My Drive\UIUC\Spring 2024\DodNDEP\Results\Data Analysis\Group A\GA [4.0] Analysis.xlsx', 'Sheet', 'MC');

% Group B Module 1-4 Data
data_PC_GB_mod_1 = readmatrix('G:\My Drive\UIUC\Spring 2024\DodNDEP\Results\Data Analysis\Group B\GB [1.0] Introduction and Rocket Hardware.xlsx', 'Sheet', 'MC');
data_PC_GB_mod_2 = readmatrix('G:\My Drive\UIUC\Spring 2024\DodNDEP\Results\Data Analysis\Group B\GB [2.0] Fundamentals of Rocketry.xlsx', 'Sheet', 'MC');
data_PC_GB_mod_3 = readmatrix('G:\My Drive\UIUC\Spring 2024\DodNDEP\Results\Data Analysis\Group B\GB [3.0] Modeling Rocket Mechanics.xlsx', 'Sheet', 'MC');
data_PC_GB_mod_4 = readmatrix('G:\My Drive\UIUC\Spring 2024\DodNDEP\Results\Data Analysis\Group B\GB [4.0] Analysis.xlsx', 'Sheet', 'MC');

num_stud_GA = length(data_PC_GA_mod_1);
num_stud_GB = length(data_PC_GB_mod_1);

% Module 1-4 Solutions
sol_mod_1 = [4, 4, 1, 2, 2, 4, 1, 1100]';
sol_mod_2 = [2, 3, 1, 3, 3, 1]';
sol_mod_3 = [2, 11010, 2, 6, 40, 100]';
sol_mod_4 = [11011, 2, 4, 1, 1]';

% Group A solution matrix for modules 1-4
sol_mat_GA_mod_1 = repmat(sol_mod_1,1,num_stud_GA)'; 
sol_mat_GA_mod_2 = repmat(sol_mod_2,1,num_stud_GA)'; 
sol_mat_GA_mod_3 = repmat(sol_mod_3,1,num_stud_GA)'; 
sol_mat_GA_mod_4 = repmat(sol_mod_4,1,num_stud_GA)'; 

% Group B solution matrix for modules 1-4
sol_mat_GB_mod_1 = repmat(sol_mod_1,1,num_stud_GB)'; 
sol_mat_GB_mod_2 = repmat(sol_mod_2,1,num_stud_GB)'; 
sol_mat_GB_mod_3 = repmat(sol_mod_3,1,num_stud_GB)'; 
sol_mat_GB_mod_4 = repmat(sol_mod_4,1,num_stud_GB)';  

% Group A difference in scores
diff_GA_mod_1 = abs(sol_mat_GA_mod_1 - data_PC_GA_mod_1);
diff_GA_mod_2 = abs(sol_mat_GA_mod_2 - data_PC_GA_mod_2);
diff_GA_mod_3 = abs(sol_mat_GA_mod_3 - data_PC_GA_mod_3);
diff_GA_mod_4 = abs(sol_mat_GA_mod_4 - data_PC_GA_mod_4);

% Group B difference in scores
diff_GB_mod_1 = abs(sol_mat_GB_mod_1 - data_PC_GB_mod_1);
diff_GB_mod_2 = abs(sol_mat_GB_mod_2 - data_PC_GB_mod_2);
diff_GB_mod_3 = abs(sol_mat_GB_mod_3 - data_PC_GB_mod_3);
diff_GB_mod_4 = abs(sol_mat_GB_mod_4 - data_PC_GB_mod_4);

% Number of students x Number of questions for Group A and B
[m_GA, n_GA] = size(sol_mat_GA_mod_1);
[m_GB, n_GB] = size(sol_mat_GB_mod_1);

% Creating vectors to store scores
score_GA_mod_1 = zeros(m_GA,1);
score_GA_mod_2 = score_GA_mod_1;
score_GA_mod_3 = score_GA_mod_1;
score_GA_mod_4 = score_GA_mod_1;

score_GB_mod_1 = zeros(m_GB,1);
score_GB_mod_2 = score_GB_mod_1;
score_GB_mod_3 = score_GB_mod_1;
score_GB_mod_4 = score_GB_mod_1;

% Calculate scores for Group A and B for module 1-4
for k = 1:m_GA
    score_GA_mod_1(k) = 100*sum(diff_GA_mod_1(k, :) == 0)/n_GA;
    score_GA_mod_2(k) = 100*sum(diff_GA_mod_2(k, :) == 0)/n_GA;
    score_GA_mod_3(k) = 100*sum(diff_GA_mod_3(k, :) == 0)/n_GA;
    score_GA_mod_4(k) = 100*sum(diff_GA_mod_4(k, :) == 0)/n_GA;
end

for j = 1:m_GB
    score_GB_mod_1(j) = 100*sum(diff_GB_mod_1(j, :) == 0)/n_GB;
    score_GB_mod_2(j) = 100*sum(diff_GB_mod_2(j, :) == 0)/n_GB;
    score_GB_mod_3(j) = 100*sum(diff_GB_mod_3(j, :) == 0)/n_GB;
    score_GB_mod_4(j) = 100*sum(diff_GB_mod_4(j, :) == 0)/n_GB;
end

x_tick_values = [0.8, 2, 3.2, 4.2];
x_label = {'Module 1', 'Module 2', 'Module 3', 'Module 4'};

figure;
pos_1 = [1, 2, 3, 4];
pos_2 = pos_1 + 0.5;

mean_GA = [score_GA_mod_1, score_GA_mod_2, score_GA_mod_3, score_GA_mod_4];
mean_GB = [score_GB_mod_1, score_GB_mod_2, score_GB_mod_3, score_GB_mod_4];

h1 = boxplot(mean_GA, 'positions', pos_1, 'colors', 'b');
hold on
h2 = boxplot(mean_GB, 'positions', pos_2, 'colors', 'r');

% Example of adjusting XTickLabel location
set(gca, 'XTickLabel', x_label, 'XTick', x_tick_values, 'FontSize', 24, 'FontName', 'Times New Roman');

% Modifying line widths
set(findobj(h1, 'type', 'line'), 'LineWidth', 1.5); 
set(findobj(h2, 'type', 'line'), 'LineWidth', 1.5); 

% Add mean values to the plot
for i = 1:numel(pos_1)
    text(pos_1(i), mean(mean_GA(:,i)) + 1, sprintf('\\mu = %.2f', mean(mean_GA(:,i))), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'b');
end
for i = 1:numel(pos_2)
    text(pos_2(i), mean(mean_GB(:,i)) + 1, sprintf('\\mu = %.2f', mean(mean_GB(:,i))), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'Color', 'r');
end

axis square; % Make axis square
grid on; % Turn on grid
title('Mid-Course Box Plot')
ylabel('Score %')
hold off
