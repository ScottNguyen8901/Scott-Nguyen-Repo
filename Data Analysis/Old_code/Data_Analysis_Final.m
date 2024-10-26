sclear, clc, close all

data_extract_GA = importdata('Group A\GA [8.0] Total.xlsx');
data_extract_GB = importdata('Group B\GB [8.0] Total.xlsx');

data_int_EC_GA = readmatrix('Group A\Interest Assessment.xlsx', 'Sheet', 'EC');
data_int_EC_GB = readmatrix('Group B\Interest Assessment.xlsx', 'Sheet', 'EC');

data_se_EC_GA = readmatrix('Group A\Self-Efficacy.xlsx', 'Sheet', 'EC');
data_se_EC_GB = readmatrix('Group B\Self-Efficacy.xlsx', 'Sheet', 'EC');

data_EC_GA = data_extract_GA.data.EC;
data_EC_GB = data_extract_GB.EC;

sol = [4, 4, 1, 2, 2, 4, 1, 1100, 2, 3, 1, 3, 3, 1, 2, 011010, 2, 6, 40, 100, 11011, 2, 4, 1, 1]';
sol_mat_GB = repmat(sol,1,length(data_EC_GB))'; 
sol_mat_GA = repmat(sol,1,length(data_EC_GA))'; 

[m_GB, n_GB] = size(sol_mat_GB);
[m_GA, n_GA] = size(sol_mat_GA);

diff_EC_GA = abs(sol_mat_GA - data_EC_GA);
diff_EC_GB = abs(sol_mat_GB - data_EC_GB);

score_EC_GA = zeros(25,1);
score_EC_GB = zeros(25,1);

for k = 1:m_GB
    score_EC_GB(k) = 100*sum(diff_EC_GB(k, :) == 0)/n_GB;
end

for k = 1:m_GA
    score_EC_GA(k) = 100*sum(diff_EC_GA(k, :) == 0)/n_GA;
end

mean_score_GA =  mean(score_EC_GA,2);
mean_score_GB = mean(score_EC_GB,2);

% Perform the t-tests
[h_EC, p_EC, ~] = ttest2(mean(data_EC_GA,2), mean(data_EC_GB,2), 'alpha', 0.01, 'Tail', 'left');
[h_int, p_int, ~] = ttest2(mean(data_int_EC_GA,2), mean(data_int_EC_GB,2), 'alpha', 0.01, 'Tail', 'right');
[h_se, p_se, ~] = ttest2(mean(data_se_EC_GA,2), mean(data_se_EC_GB,2), 'alpha', 0.01, 'Tail', 'left');

% Create a table to store results
table_tests = cell(4,3); % Initialize the table with appropriate size

% Headers
table_tests{1,1} = ''; % Empty cell for formatting
table_tests{1,2} = 'P-Value';
table_tests{1,3} = 'Decision';

% Labels for each row
table_tests{2,1} = 'Final';
table_tests{3,1} = 'Interest';
table_tests{4,1} = 'Self-Efficacy';

% Fill in the table with data
table_tests{2,2} = p_EC;
table_tests{3,2} = p_int;
table_tests{4,2} = p_se;

% Decision based on p-value < alpha level (0.01)
alpha = 0.01;
table_tests{2,3} = decision(p_EC, alpha);
table_tests{3,3} = decision(p_int, alpha);
table_tests{4,3} = decision(p_se, alpha);

% Print the table
fprintf('\nHypothesis Test Results:\n');
fprintf('%-15s %-15s %-15s\n', table_tests{1,:});
for i = 2:size(table_tests,1)
    fprintf('%-15s %-15.6e %-15s\n', table_tests{i,1}, table_tests{i,2}, table_tests{i,3});
end
%% Plotting

fs = 24; % Font size for labels and text
figure;

% Define the positions for the subplots
left = 0.1; % Left margin
bottom = 0.1; % Bottom margin
width = 0.35; % Width of each subplot
height = 0.8; % Height of each subplot
space = 0.1; % Space between subplots

% First subplot for the original data
subplot('Position', [left bottom width height]);

% Define x-coordinates for the boxplots
x1 = 1; % Position for mean_score_GA
x2 = 1.15; % Position for mean_score_GB

h1 = boxplot(mean_score_GA, 'positions', x1, 'colors', 'b');
hold on
h2 = boxplot(mean_score_GB, 'positions', x2, 'colors', 'r');

% Modifying line widths
set(findobj(h1, 'type', 'line'), 'LineWidth', 1.5); 
set(findobj(h2, 'type', 'line'), 'LineWidth', 1.5); 

text(x1, mean(mean_score_GA) + 5, sprintf('\\mu = %.2f', mean(mean_score_GA)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 16, 'Color', 'b', 'FontName', 'Times New Roman');
text(x2, mean(mean_score_GB) - 5, sprintf('\\mu = %.2f', mean(mean_score_GB)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 16, 'Color', 'r', 'FontName', 'Times New Roman');
ylabel('Score %', 'FontSize', fs, 'FontName', 'Times New Roman');
xlabel('')

% Add 'Final Assessment' label
set(gca, 'XTick', [mean([x1, x2])], 'XTickLabel', {'Final Assessment'}, 'FontSize', fs, 'FontName', 'Times New Roman');

axis square; % Make axis square
grid on; % Turn on grid

hold off

% Assuming mean_se_GA and mean_se_GB are the mean values of self-efficacy for Group A and Group B
% You would need to calculate these or replace them with the appropriate variables
mean_int_GA = mean(data_int_EC_GA, 2, 'omitnan');
mean_int_GB = mean(data_int_EC_GB, 2, 'omitnan');
mean_se_GA = mean(data_se_EC_GA, 2, 'omitnan');
mean_se_GB = mean(data_se_EC_GB, 2, 'omitnan');

% Define positions for the boxplots
pos_1 = 1; % Position for mean_int_GA
pos_2 = 1.5; % Position for mean_int_GB
pos_3 = 2; % Position for mean_se_GA
pos_4 = 2.5; % Position for mean_se_GB

% Define x-tick labels and their positions
x_label = {'Interest', 'Self-Efficacy'};
x_tick_values = [mean([pos_1, pos_2]), mean([pos_3, pos_4])]; % Midpoint between positions

% Second subplot for the interest and self-efficacy data
subplot('Position', [left+width+space bottom width height]);

% Boxplot for interest assessment data for Group A
h3 = boxplot(mean_int_GA, 'positions', pos_1, 'colors', 'b', 'Widths', 0.5); % Increase the width here
hold on
% Boxplot for interest assessment data for Group B
h4 = boxplot(mean_int_GB, 'positions', pos_2, 'colors', 'r', 'Widths', 0.5); % And here

% Boxplot for self-efficacy data for Group A
h5 = boxplot(mean_se_GA, 'positions', pos_3, 'colors', 'b', 'Widths', 0.5); % And here
% Boxplot for self-efficacy data for Group B
h6 = boxplot(mean_se_GB, 'positions', pos_4, 'colors', 'r', 'Widths', 0.5); % And here

% Modifying line widths
set(findobj(gca, 'type', 'line'), 'LineWidth', 1.5); 

% Add mean values to the plot for interest
text(pos_1, mean(mean_int_GA) + 0.03 * range(ylim), sprintf('\\mu = %.2f', mean(mean_int_GA)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 15, 'Color', 'b');
text(pos_2, mean(mean_int_GB) + 0.03 * range(ylim), sprintf('\\mu = %.2f', mean(mean_int_GB)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 15, 'Color', 'r');

% Add mean values to the plot for self-efficacy
text(pos_3, mean(mean_se_GA) + 0.03 * range(ylim), sprintf('\\mu = %.2f', mean(mean_se_GA)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 15, 'Color', 'b');
text(pos_4, mean(mean_se_GB) - 0.01 * range(ylim), sprintf('\\mu = %.2f', mean(mean_se_GB)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 15, 'Color', 'r');


% Example of adjusting XTickLabel location
set(gca, 'XTickLabel', x_label, 'XTick', x_tick_values, 'FontSize', 24, 'FontName', 'Times New Roman');

ylabel('Score', 'FontSize', 24, 'FontName', 'Times New Roman')
axis square; % Make axis square
grid on; % Turn on grid
legend([h3(1), h4(1)], 'Group A', 'Group B', 'TextColor', 'black', 'Location', 'best');

hold off

