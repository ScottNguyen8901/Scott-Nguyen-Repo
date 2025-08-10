clear, clc, close all

% Group A data
data_PC_GA = readmatrix('G:\My Drive\UIUC\Spring 2024\DodNDEP\Results\Data Analysis\Group A\GA [3.0] Modeling Rocket Mechanics.xlsx', 'Sheet', 'PC');
data_MC_GA = readmatrix('G:\My Drive\UIUC\Spring 2024\DodNDEP\Results\Data Analysis\Group A\GA [3.0] Modeling Rocket Mechanics.xlsx', 'Sheet', 'MC');
data_EC_GA = readmatrix('G:\My Drive\UIUC\Spring 2024\DodNDEP\Results\Data Analysis\Group A\GA [3.0] Modeling Rocket Mechanics.xlsx', 'Sheet', 'EC');

% Group B data
data_PC_GB = readmatrix('G:\My Drive\UIUC\Spring 2024\DodNDEP\Results\Data Analysis\Group B\GB [3.0] Modeling Rocket Mechanics.xlsx', 'Sheet', 'PC');
data_MC_GB = readmatrix('G:\My Drive\UIUC\Spring 2024\DodNDEP\Results\Data Analysis\Group B\GB [3.0] Modeling Rocket Mechanics.xlsx', 'Sheet', 'MC');
data_EC_GB = readmatrix('G:\My Drive\UIUC\Spring 2024\DodNDEP\Results\Data Analysis\Group B\GB [3.0] Modeling Rocket Mechanics.xlsx', 'Sheet', 'EC');

sol = [2, 011010, 2, 6, 40, 100]';
sol_mat_GB = repmat(sol,1,length(data_PC_GB))'; 
sol_mat_GA = repmat(sol,1,length(data_PC_GA))'; 

[m_GB, n_GB] = size(sol_mat_GB);
[m_GA, n_GA] = size(sol_mat_GA);

diff_PC_GB = abs(sol_mat_GB - data_PC_GB);
diff_MC_GB = abs(sol_mat_GB - data_MC_GB);
diff_EC_GB = abs(sol_mat_GB - data_EC_GB);

diff_PC_GA = abs(sol_mat_GA - data_PC_GA);
diff_MC_GA = abs(sol_mat_GA - data_MC_GA);
diff_EC_GA = abs(sol_mat_GA - data_EC_GA);

score_PC_GB = zeros(25,1);
score_MC_GB = zeros(25,1);
score_EC_GB = zeros(25,1);

score_PC_GA = zeros(25,1);
score_MC_GA = zeros(25,1);
score_EC_GA = zeros(25,1);

for k = 1:m_GB
    score_PC_GB(k) = 100*sum(diff_PC_GB(k, :) == 0)/n_GB;
    score_MC_GB(k) = 100*sum(diff_MC_GB(k, :) == 0)/n_GB;
    score_EC_GB(k) = 100*sum(diff_EC_GB(k, :) == 0)/n_GB;
end

for k = 1:m_GA
    score_PC_GA(k) = 100*sum(diff_PC_GA(k, :) == 0)/n_GA;
    score_MC_GA(k) = 100*sum(diff_MC_GA(k, :) == 0)/n_GA;
    score_EC_GA(k) = 100*sum(diff_EC_GA(k, :) == 0)/n_GA;
end

MC_to_PC_GB = score_MC_GB - score_PC_GB;
EC_to_MC_GB = score_EC_GB - score_MC_GB;
EC_to_PC_GB = score_EC_GB - score_PC_GB;

MC_to_PC_GA = score_MC_GA - score_PC_GA;
EC_to_MC_GA = score_EC_GA - score_MC_GA;
EC_to_PC_GA = score_EC_GA - score_PC_GA;

med_PC_GA = median(score_PC_GA,"all");
med_MC_GA = median(score_MC_GA,"all");
med_EC_GA = median(score_EC_GA,"all");

med_PC_GB = median(score_PC_GB,"all");
med_MC_GB = median(score_MC_GB,"all");
med_EC_GB = median(score_EC_GB,"all");

mean_PC_GB = mean(score_PC_GB,"all");
mean_MC_GB = mean(score_MC_GB,"all");
mean_EC_GB = mean(score_EC_GB,"all");

mean_PC_GA = mean(score_PC_GA,"all");
mean_MC_GA = mean(score_MC_GA,"all");
mean_EC_GA = mean(score_EC_GA,"all");

mean_val_GA  = [mean_PC_GA, mean_MC_GA , mean_EC_GA];
mean_val_GB  = [mean_PC_GB, mean_MC_GB , mean_EC_GB];

mean_score_GA = [mean(score_PC_GA,2), mean(score_MC_GA,2), mean(score_EC_GA,2)];
mean_score_GB = [mean(score_PC_GB,2), mean(score_MC_GB,2), mean(score_EC_GB,2)];

med_val_GA = [med_PC_GA, med_MC_GA, med_EC_GA];
med_val_GB = [med_PC_GB, med_MC_GB, med_EC_GB];

mean_MC_to_PC_GA = mean_MC_GA - mean_PC_GA;
mean_EC_to_MC_GA = mean_EC_GA - mean_MC_GA;
mean_EC_to_PC_GA = mean_EC_GA - mean_PC_GA;

mean_MC_to_PC_GB = mean_MC_GB - mean_PC_GB;
mean_EC_to_MC_GB = mean_EC_GB - mean_MC_GB;
mean_EC_to_PC_GB = mean_EC_GB - mean_PC_GB;

med_MC_to_PC_GA = med_MC_GA - med_PC_GA;
med_EC_to_MC_GA = med_EC_GA - med_MC_GA;
med_EC_to_PC_GA = med_EC_GA - med_PC_GA;

med_MC_to_PC_GB = med_MC_GB - med_PC_GB;
med_EC_to_MC_GB = med_EC_GB - med_MC_GB;
med_EC_to_PC_GB = med_EC_GB - med_PC_GB;

mean_val_diff_GA = [mean_MC_to_PC_GA, mean_EC_to_MC_GA, mean_EC_to_PC_GA];
mean_val_diff_GB = [mean_MC_to_PC_GB, mean_EC_to_MC_GB, mean_EC_to_PC_GB];

med_val_diff_GA = [med_MC_to_PC_GA, med_EC_to_MC_GA, med_EC_to_PC_GA];
med_val_diff_GB = [med_MC_to_PC_GB, med_EC_to_MC_GB, med_EC_to_PC_GB];

scale = {'Strongly Disagree [1]', 'Disagree[2]', 'Neutral[3]', 'Agree[4]', 'Strongly Agree[5]'};
xlab = {'Pre-Course -> Mid-Course', 'Mid-Course -> End-Course', 'Pre-Course -> End-Course'};

[H_PC_GA, ~, ~] = swtest(score_PC_GA, 0.01);
[H_MC_GA, ~, ~] = swtest(score_MC_GA, 0.01);
[H_EC_GA, ~, ~] = swtest(score_EC_GA, 0.01);

[H_PC_GB, ~, ~] = swtest(score_PC_GB, 0.01);
[H_MC_GB, ~, ~] = swtest(score_MC_GB, 0.01);
[H_EC_GB, ~, ~] = swtest(score_EC_GB, 0.01);

if H_PC_GA == 0 && H_PC_GB == 0
    % Perform ttest2
    [h_PC, p_PC, stats_PC] = ttest2(mean(data_PC_GA,2), mean(data_PC_GB,2), 'alpha', 0.01);
    test_PC = 0;
else
    % Perform ranksum test
    [p_PC, h_PC, stats_PC] = ranksum(median(data_PC_GA,2), median(data_PC_GB,2), 'tail', 'both', 'alpha', 0.01);
    test_PC = 1;
end

if H_MC_GA == 0 && H_MC_GB == 0
    % Perform ttest2
    [h_MC, p_MC, stats_MC] = ttest2(mean(data_MC_GA,2), mean(data_MC_GB,2), 'alpha', 0.01);
    test_MC = 0;
else
    % Perform ranksum test
    [p_MC, h_MC, stats_MC] = ranksum(median(data_MC_GA,2), median(data_MC_GB,2), 'tail', 'both', 'alpha', 0.01);
    test_MC = 1;
end

if H_EC_GA == 0 && H_EC_GB == 0
    % Perform ttest2
    [h_EC, p_EC, stats_EC] = ttest2(mean(data_EC_GA,2), mean(data_EC_GB,2), 'alpha', 0.01);
    test_EC = 0;
else
    % Perform ranksum test
    [p_EC, h_EC, stats_EC] = ranksum(median(data_EC_GA,2), median(data_EC_GB,2), 'tail', 'left', 'alpha', 0.01);
    test_EC = 1;
end

[H_PCMC_GA, p_PCMC_GA, ~] = ttest2(mean(data_PC_GA,2), mean(data_MC_GA,2), 'tail', 'left', 'alpha', 0.01);
[H_PCMC_GB, p_PCMC_GB, ~] = ttest2(mean(data_PC_GB,2), mean(data_MC_GB,2), 'tail', 'left', 'alpha', 0.01);

[H_MCEC_GA, p_MCEC_GA, ~] = ttest2(mean(data_MC_GA,2), mean(data_EC_GA,2), 'tail', 'left', 'alpha', 0.01);
[H_MCEC_GB, p_MCEC_GB, ~] = ttest2(mean(data_MC_GB,2), mean(data_EC_GB,2), 'tail', 'left', 'alpha', 0.01);

[H_PCEC_GA, p_PCEC_GA, ~] = ttest2(mean(data_PC_GA,2), mean(data_EC_GA,2), 'tail', 'left', 'alpha', 0.01);
[H_PCEC_GB, p_PCEC_GB, ~] = ttest2(mean(data_PC_GB,2), mean(data_EC_GB,2), 'tail', 'left', 'alpha', 0.01);

data_analysis = {'P-Value', 'Alpha Value', 'Test';
                 p_PC, '0.01', ternary(test_PC, 'T-test', 'Wilcoxon');
                 p_MC, '0.01', ternary(test_MC, 'T-test', 'Wilcoxon');
                 p_EC, '0.01', ternary(test_EC, 'T-test', 'Wilcoxon')};

% Create a table to store results
results_table = {'', 'P-Value', 'Alpha Level', 'Test', 'Reject Null'};
results_table{2,1} = 'PC';
results_table{3,1} = 'MC';
results_table{4,1} = 'EC';

% Fill in the table with data
results_table{2,2} = data_analysis{2,1};
results_table{3,2} = data_analysis{3,1};
results_table{4,2} = data_analysis{4,1};

results_table{2,3} = data_analysis{2,2};
results_table{3,3} = data_analysis{3,2};
results_table{4,3} = data_analysis{4,2};

results_table{2,4} = data_analysis{2,3};
results_table{3,4} = data_analysis{3,3};
results_table{4,4} = data_analysis{4,3};

results_table{2,5} = ternary(h_PC, 'Yes', 'No');
results_table{3,5} = ternary(h_MC, 'Yes', 'No');
results_table{4,5} = ternary(h_EC, 'Yes', 'No');

% Print the table
% Print the table
fprintf('Hypothesis Test Results:\n');
fprintf('%-10s %-15s %-15s %-15s %-15s\n', results_table{1,:});
for i = 2:size(results_table,1)
    fprintf('%-10s %-15.6e %-15.2f %-15s %-15s\n', results_table{i,1}, results_table{i,2}, str2double(results_table{i,3}), results_table{i,4}, results_table{i,5});
end

results_table_comp = {'Group', 'Mean PC', 'Mean MC', 'Mean EC', 'Mean PC->MC', 'Mean MC->EC', 'Mean PC->EC'};

results_table_comp{2,1} = 'Group A';
results_table_comp{3,1} = 'Group B';

results_table_comp{2,2} = mean_PC_GA;
results_table_comp{3,2} = mean_PC_GB;

results_table_comp{2,3} = mean_MC_GA;
results_table_comp{3,3} = mean_MC_GB;

results_table_comp{2,4} = mean_EC_GA;
results_table_comp{3,4} = mean_EC_GB;

results_table_comp{2,5} = mean_MC_to_PC_GA;
results_table_comp{3,5} = mean_MC_to_PC_GB;

results_table_comp{2,6} = mean_EC_to_MC_GA;
results_table_comp{3,6} = mean_EC_to_MC_GB;

results_table_comp{2,7} = mean_EC_to_PC_GA;
results_table_comp{3,7} = mean_EC_to_PC_GB;

% Print the table
fprintf('\nComparison Mean Results:\n');
fprintf('%-12s %-12s %-12s %-12s %-12s %-12s %-12s\n', results_table_comp{1,:});
for i = 2:size(results_table_comp,1)
    fprintf('%-12s %-12.2f %-12.2f %-12.2f %-12.2f %-12.2f %-12.2f\n', results_table_comp{i,:});
end

results_table_comp_med = {'Group', 'Median PC', 'Median MC', 'Median EC', 'Median PC->MC', 'Median MC->EC', 'Median PC->EC'};

results_table_comp_med{2,1} = 'Group A';
results_table_comp_med{3,1} = 'Group B';

results_table_comp_med{2,2} = med_PC_GA;
results_table_comp_med{3,2} = med_PC_GB;

results_table_comp_med{2,3} = med_MC_GA;
results_table_comp_med{3,3} = med_MC_GB;

results_table_comp_med{2,4} = med_EC_GA;
results_table_comp_med{3,4} = med_EC_GB;

results_table_comp_med{2,5} = med_MC_to_PC_GA;
results_table_comp_med{3,5} = med_MC_to_PC_GB;

results_table_comp_med{2,6} = med_EC_to_MC_GA;
results_table_comp_med{3,6} = med_EC_to_MC_GB;

results_table_comp_med{2,7} = med_EC_to_PC_GA;
results_table_comp_med{3,7} = med_EC_to_PC_GB;

% Print the table
fprintf('\nComparison Median Results:\n');
fprintf('%-12s %-12s %-12s %-12s %-12s %-12s %-12s\n', results_table_comp_med{1,:});
for i = 2:size(results_table_comp,1)
    fprintf('%-12s %-12.2f %-12.2f %-12.2f %-12.2f %-12.2f %-12.2f\n', results_table_comp_med{i,:});
end

% Perform hypothesis tests and store results for Group A
test_results_A = {'P-Value', 'Alpha Value', 'Test', 'Reject Null';
                  p_PCMC_GA, '0.01', 'T-test', ternary(H_PCMC_GA, 'Yes', 'No');
                  p_MCEC_GA, '0.01', 'T-test', ternary(H_MCEC_GA, 'Yes', 'No');
                  p_PCEC_GA, '0.01', 'T-test', ternary(H_PCEC_GA, 'Yes', 'No')};

% Create a table to store results for Group A
table_A = cell(4,5); % Initialize the table with appropriate size
table_A{1,1} = ''; % Empty cell for formatting

% Headers for Group A
table_A{1,2} = 'P-Value';
table_A{1,3} = 'Alpha Level';
table_A{1,4} = 'Test';
table_A{1,5} = 'Reject Null';

% Labels for each row for Group A
table_A{2,1} = 'PC->MC';
table_A{3,1} = 'MC->EC';
table_A{4,1} = 'PC->EC';

% Fill in the table with data for Group A
for i = 2:size(table_A,1)
    table_A{i,2} = test_results_A{i,1};
    table_A{i,3} = test_results_A{i,2};
    table_A{i,4} = test_results_A{i,3};
    table_A{i,5} = test_results_A{i,4};
end

% Print the table for Group A
fprintf('Hypothesis Test Results for Group A:\n');
fprintf('%-10s %-15s %-15s %-15s %-15s\n', table_A{1,:});
for i = 2:size(table_A,1)
    fprintf('%-10s %-15.6e %-15.2f %-15s %-15s\n', table_A{i,1}, table_A{i,2}, str2double(table_A{i,3}), table_A{i,4}, table_A{i,5});
end

% Perform hypothesis tests and store results for Group B
test_results_B = {'P-Value', 'Alpha Value', 'Test', 'Reject Null';
                  p_PCMC_GB, '0.01', 'T-test', ternary(H_PCMC_GB, 'Yes', 'No');
                  p_MCEC_GB, '0.01', 'T-test', ternary(H_MCEC_GB, 'Yes', 'No');
                  p_PCEC_GB, '0.01', 'T-test', ternary(H_PCEC_GB, 'Yes', 'No')};

% Create a table to store results for Group B
table_B = cell(4,5); % Initialize the table with appropriate size
table_B{1,1} = ''; % Empty cell for formatting

% Headers for Group B
table_B{1,2} = 'P-Value';
table_B{1,3} = 'Alpha Level';
table_B{1,4} = 'Test';
table_B{1,5} = 'Reject Null';

% Labels for each row for Group B
table_B{2,1} = 'PC->MC';
table_B{3,1} = 'MC->EC';
table_B{4,1} = 'PC->EC';

% Fill in the table with data for Group B
for i = 2:size(table_B,1)
    table_B{i,2} = test_results_B{i,1};
    table_B{i,3} = test_results_B{i,2};
    table_B{i,4} = test_results_B{i,3};
    table_B{i,5} = test_results_B{i,4};
end

% Print the table for Group B
fprintf('\nHypothesis Test Results for Group B:\n');
fprintf('%-10s %-15s %-15s %-15s %-15s\n', table_B{1,:});
for i = 2:size(table_B,1)
    fprintf('%-10s %-15.6e %-15.2f %-15s %-15s\n', table_B{i,1}, table_B{i,2}, str2double(table_B{i,3}), table_B{i,4}, table_B{i,5});
end

time = [1,2,3];
bin_edges = 0:10:100;
x_label = {'Pre-Course', 'Mid-Course', 'End-Course'};
x_label_comp = {'PC -> MC', 'MC -> EC', 'PC -> EC'};
x_tick_values = [0.8, 2, 3.2];

fs = 24;
figure;

subplot(1,3,1);
hold on;
histogram(mean_score_GA(:,1), bin_edges, 'EdgeColor', 'b', 'FaceColor', [0.5 0.5 1], 'FaceAlpha', 0.5); % GA histogram with lighter blue fill
histogram(mean_score_GB(:,1), bin_edges, 'EdgeColor', 'r', 'FaceColor', [1 0 0], 'FaceAlpha', 0.5); % GB histogram with dark red fill
hold off;
grid on
axis square
lgd = legend('', '');
title(lgd, 'Pre-Course', 'FontSize', fs, 'FontName', 'Times New Roman');
legend('Location', 'north');
ylabel('Count', 'FontSize', fs, 'FontName', 'Times New Roman');

subplot(1,3,2);
hold on;
histogram(mean_score_GA(:,2), bin_edges, 'EdgeColor', 'b', 'FaceColor', [0.5 0.5 1], 'FaceAlpha', 0.5); % GA histogram with lighter blue fill
histogram(mean_score_GB(:,2), bin_edges, 'EdgeColor', 'r', 'FaceColor', [1 0 0], 'FaceAlpha', 0.5); % GB histogram with dark red fill
hold off;
grid on
axis square
lgd = legend('', '');
title(lgd, 'Mid-Course', 'FontSize', fs, 'FontName', 'Times New Roman');
legend('Location', 'north');
xlabel('Score', 'FontSize', fs, 'FontName', 'Times New Roman');

subplot(1,3,3);
hold on;
histogram(mean_score_GA(:,3), bin_edges, 'EdgeColor', 'b', 'FaceColor', [0.5 0.5 1], 'FaceAlpha', 0.5); % GA histogram with lighter blue fill
histogram(mean_score_GB(:,3), bin_edges, 'EdgeColor', 'r', 'FaceColor', [1 0 0], 'FaceAlpha', 0.5); % GB histogram with dark red fill
hold off;
grid on
axis square
lgd = legend('GA', 'GB');
legend('Location', 'north');
title(lgd, 'End-Course', 'FontSize', fs, 'FontName', 'Times New Roman');

figure;

subplot(1,2,1)
pos_1 = [1, 2, 3] - 0.25;
pos_2 = pos_1 + 0.45;

h1 = boxplot(mean_score_GA, 'positions', pos_1, 'colors', 'b');
hold on
h2 = boxplot(mean_score_GB, 'positions', pos_2, 'colors', 'r');

% Example of adjusting XTickLabel location
set(gca, 'XTickLabel', x_label, 'XTick', x_tick_values, 'FontSize', 24, 'FontName', 'Times New Roman');

% Modifying line widths
set(findobj(h1, 'type', 'line'), 'LineWidth', 1.5); 
set(findobj(h2, 'type', 'line'), 'LineWidth', 1.5); 

% Add mean values to the plot
for i = 1:numel(pos_1)
    text(pos_1(i), mean(mean_score_GA(:,i)) +5, sprintf('\\mu = %.2f', mean(mean_score_GA(:,i))), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10, 'Color', 'b');
end
for i = 1:numel(pos_2)
    text(pos_2(i), mean(mean_score_GB(:,i)) - 7, sprintf('\\mu = %.2f', mean(mean_score_GB(:,i))), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10, 'Color', 'r');
end

ylabel('Score %')
axis square; % Make axis square
grid on; % Turn on grid

hold off

subplot(1,2,2)
alpha_level = 0.01; % Alpha level for hypothesis testing
bar(time-0.2, mean_val_diff_GA, 0.4, 'b') % Shift the first set of bars slightly to the left
hold on
bar(time+0.2, mean_val_diff_GB, 0.4, 'r') % Shift the second set of bars slightly to the right
grid on
axis square
xticks(1:3);
set(gca, 'XTickLabel', x_label_comp, 'XTick', x_tick_values, 'FontSize', 24, 'FontName', 'Times New Roman');
ylabel('Change in Score %', 'FontSize', 24, 'FontName', 'Times New Roman')
ylim([-10, 30]) % Adjust the y-axis limit as needed
legend('Group A', 'Group B', 'FontSize', 24) 

sgtitle('Modeling Rocket Dynamics Comparison', 'FontSize', 24, 'FontName', 'Times New Roman')

% Add p-values and rejection decision with text moved up and dark blue color
text(time(1), max(mean_val_diff_GA(1), mean_val_diff_GB(1)) + 3, ...
    sprintf('p = %.1e\nReject H0: %s', p_PCMC_GA, decision(p_PCMC_GA, alpha_level)), ...
    'HorizontalAlignment', 'center', 'FontSize', 14, 'Color', 'b');

text(time(2), max(mean_val_diff_GA(2), mean_val_diff_GB(2)) + 3, ...
    sprintf('p = %.1e\nReject H0: %s', p_MCEC_GA, decision(p_MCEC_GA, alpha_level)), ...
    'HorizontalAlignment', 'center', 'FontSize', 14, 'Color', 'b');

text(time(3), max(mean_val_diff_GA(3), mean_val_diff_GB(3)) + 3, ...
    sprintf('p = %.1e\nReject H0: %s', p_PCEC_GA, decision(p_PCEC_GA, alpha_level)), ...
    'HorizontalAlignment', 'center', 'FontSize', 14, 'Color', 'b');

% Add p-values and rejection decision with text moved up and dark blue color
text(time(1), -7, ...
    sprintf('p = %.1e\nReject H0: %s', p_PCMC_GB, decision(p_PCMC_GB, alpha_level)), ...
    'HorizontalAlignment', 'center', 'FontSize', 14, 'Color', 'r');

text(time(2), -7, ...
    sprintf('p = %.1e\nReject H0: %s', p_MCEC_GB, decision(p_MCEC_GB, alpha_level)), ...
    'HorizontalAlignment', 'center', 'FontSize', 14, 'Color', 'r');

text(time(3), -7, ...
    sprintf('p = %.1e\nReject H0: %s', p_PCEC_GB, decision(p_PCEC_GB, alpha_level)), ...
    'HorizontalAlignment', 'center', 'FontSize', 14, 'Color', 'r');
    
% Add bar values on top of the bars with the Greek letter delta
for i = 1:length(time)
    text(time(i), -3, sprintf('\\Delta=%.2f', mean_val_diff_GA(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14, 'Color', 'b', 'Interpreter', 'tex');
    text(time(i), -5, sprintf('\\Delta=%.2f', mean_val_diff_GB(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 14, 'Color', 'r', 'Interpreter', 'tex');
end

