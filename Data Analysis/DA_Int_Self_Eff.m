clear, clc, close all

data_int_PC_GA = readmatrix('Group A\Interest Assessment.xlsx', 'Sheet', 'PC');
data_int_PC_GB = readmatrix('Group B\Interest Assessment.xlsx', 'Sheet', 'PC');

data_int_MC_GA = readmatrix('Group A\Interest Assessment.xlsx', 'Sheet', 'MC');
data_int_MC_GB = readmatrix('Group B\Interest Assessment.xlsx', 'Sheet', 'MC');

data_se_PC_GA = readmatrix('Group A\Self-Efficacy.xlsx', 'Sheet', 'PC');
data_se_PC_GB = readmatrix('Group B\Self-Efficacy.xlsx', 'Sheet', 'PC');

data_se_MC_GA = readmatrix('Group A\Self-Efficacy.xlsx', 'Sheet', 'MC');
data_se_MC_GB = readmatrix('Group B\Self-Efficacy.xlsx', 'Sheet', 'MC');

mean_PC_GA = [mean(data_int_PC_GA,2), mean(data_se_PC_GA,2)];
mean_PC_GB = [mean(data_int_PC_GB,2) , mean(data_se_PC_GB,2)];

mean_MC_GA = [mean(data_int_MC_GA,2), mean(data_se_MC_GA,2)];
mean_MC_GB = [mean(data_int_MC_GB,2) , mean(data_se_MC_GB,2)];

change_PC_MC_GA = mean_MC_GA - mean_PC_GA;
change_PC_MC_GB = mean_MC_GB - mean_PC_GB;

[h_int_PC_MC_GA, p_int_PC_MC_GA, ~] = ttest2(mean_PC_GA(:, 1), mean_MC_GA(:, 1), 'alpha', 0.01, 'tail', 'right');
[h_se_PC_MC_GA, p_se_PC_MC_GA, ~] = ttest2(mean_PC_GA(:, 2), mean_MC_GA(:, 2), 'alpha', 0.01, 'tail', 'left');

[h_int_PC_MC_GB, p_int_PC_MC_GB, ~] = ttest2(mean_PC_GB(:, 1), mean_MC_GB(:, 1), 'alpha', 0.01, 'tail', 'left');
[h_se_PC_MC_GB, p_se_PC_MC_GB, ~] = ttest2(mean_PC_GB(:, 2), mean_MC_GB(:, 2), 'alpha', 0.01, 'tail', 'left');

[h_int_PC, p_int_PC, ~] = ttest2(mean_PC_GA(:, 1), mean_PC_GB(:, 1), 'alpha', 0.01, 'tail', 'right');
[h_int_MC, p_int_MC, ~] = ttest2(mean_MC_GA(:, 1), mean_MC_GB(:, 1), 'alpha', 0.01);

[h_se_PC, p_se_PC, ~] = ttest2(mean_PC_GA(:, 2), mean_PC_GB(:, 2), 'alpha', 0.01, 'tail', 'left');
[h_se_MC, p_se_MC, ~] = ttest2(mean_MC_GA(:, 2), mean_MC_GB(:, 2), 'alpha', 0.01);

mean_int_GA = [mean_PC_GA(:,1), mean_MC_GA(:,1)];
mean_int_GB = [mean_PC_GB(:,1), mean_MC_GB(:,1)];

mean_se_GA = [mean_PC_GA(:,2), mean_MC_GA(:,2)];
mean_se_GB = [mean_PC_GB(:,2), mean_MC_GB(:,2)];

x_label = {'Pre-Content', 'Mid-Content'};
x_tick_values = [1, 2];
pos_1 = [1, 2] - 0.25;
pos_2 = pos_1 + 0.45;

figure;
subplot(1,2,1)

subplot(1,2,1)
h1 = boxplot(mean_int_GA, 'positions', pos_1, 'colors', 'b', 'BoxWidth', 0.6); % Bigger box
hold on
h2 = boxplot(mean_int_GB, 'positions', pos_2, 'colors', 'r', 'BoxWidth', 0.6); % Bigger box

outliers_h1 = findobj(h1,'tag','Outliers');
set(outliers_h1, 'MarkerEdgeColor', 'b');

% Example of adjusting XTickLabel location
set(gca, 'XTickLabel', x_label, 'XTick', x_tick_values, 'FontSize', 24, 'FontName', 'Times New Roman');

% Modifying line widths
set(findobj(h1, 'type', 'line'), 'LineWidth', 1.5); 
set(findobj(h2, 'type', 'line'), 'LineWidth', 1.5); 

% Add mean values to the plot
for i = 1:numel(pos_1)
    % Calculate the mean ignoring NaN values
    mu_GA = mean(mean_int_GA(:,i));
    % Display the mean value on the plot
    text(pos_1(i), mu_GA + 0.03 * range(ylim), sprintf('\\mu = %.2f', mu_GA), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 15, 'Color', 'b');
end

for i = 1:numel(pos_2)
    % Calculate the mean ignoring NaN values
    mu_GB = mean(mean_int_GB(:,i));
    % Display the mean value on the plot
    text(pos_2(i), mu_GB + 0.03 * range(ylim), sprintf('\\mu = %.2f', mu_GB), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 15, 'Color', 'r');
end

ylabel('Score', 'FontSize', 24, 'FontName', 'Times New Roman')
axis square; % Make axis square
grid on; % Turn on grid

hold off

subplot(1,2,2)
h3 = boxplot(mean_se_GA, 'positions', pos_1, 'colors', 'b');
hold on
h4 = boxplot(mean_se_GB, 'positions', pos_2, 'colors', 'r');

outliers_h3 = findobj(h3,'tag','Outliers');
set(outliers_h3, 'MarkerEdgeColor', 'b');

% Example of adjusting XTickLabel location
set(gca, 'XTickLabel', x_label, 'XTick', x_tick_values, 'FontSize', 24, 'FontName', 'Times New Roman');

% Modifying line widths
set(findobj(h3, 'type', 'line'), 'LineWidth', 1.5); 
set(findobj(h4, 'type', 'line'), 'LineWidth', 1.5); 

% Add mean values to the plot
for i = 1:numel(pos_1)
    % Calculate the mean ignoring NaN values
    mu_GA = mean(mean_se_GA(:,i));
    % Display the mean value on the plot
    text(pos_1(i), mu_GA + 0.03 * range(ylim), sprintf('\\mu = %.2f', mu_GA), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 15, 'Color', 'b');
end

for i = 1:numel(pos_2)
    % Calculate the mean ignoring NaN values
    mu_GB = mean(mean_se_GB(:,i));
    % Display the mean value on the plot
    text(pos_2(i), mu_GB + 0.02 * range(ylim), sprintf('\\mu = %.2f', mu_GB), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 15, 'Color', 'r');
end

ylabel('Score', 'FontSize', 24, 'FontName', 'Times New Roman')
axis square; % Make axis square
grid on; % Turn on grid
legend([h1(1), h2(1)], 'Group A', 'Group B', 'TextColor', 'black', 'Location', 'best');

hold off

% Define the font size variable
fs = 24; 
alpha_level = 0.01; % Alpha level for hypothesis testing

% Define the time points for the x-axis
time = [1, 2];

% Define the labels for the x-axis
x_label_comp = {'Interest', 'Self-Efficacy'};
x_tick_values = [1, 2];

% Calculate the mean value differences for Group A and Group B
mean_val_diff_GA = mean(change_PC_MC_GA, 1); % Mean across rows for each column
mean_val_diff_GB = mean(change_PC_MC_GB, 1); % Mean across rows for each column

% Plot the bar graph for the mean value differences
figure;
bar(time-0.2, mean_val_diff_GA, 0.4, 'b') % Shift the first set of bars slightly to the left
hold on
bar(time+0.2, mean_val_diff_GB, 0.4, 'r') % Shift the second set of bars slightly to the right
grid on
axis square
xticks(time);
set(gca, 'XTickLabel', x_label_comp, 'XTick', x_tick_values, 'FontSize', 24, 'FontName', 'Times New Roman');
ylim([-1.2, 1]) % Adjust the y-axis limit as needed
ylabel('Change in Score', 'FontSize', 24, 'FontName', 'Times New Roman')
legend('Group A', 'Group B', 'FontSize', 24) 

% Add the hypothesis testing results to the plot
% Since p_int_PC_MC_GA and p_int_PC_MC_GB are scalars, we do not use a loop or index them

% Add mean value differences on top of the bars for Group A
text(time(1) - 0.2, + 0.05, sprintf('\\Delta=%.2f', mean_val_diff_GA(1)), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fs, 'Color', 'b', 'FontName', 'Times New Roman', 'Interpreter', 'tex');
text(time(2) - 0.2, mean_val_diff_GA(2) + 0.05, sprintf('\\Delta=%.2f', mean_val_diff_GA(2)), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fs, 'Color', 'b', 'FontName', 'Times New Roman', 'Interpreter', 'tex');

% Add mean value differences on top of the bars for Group B
text(time(1) + 0.2, mean_val_diff_GB(1) + 0.05, sprintf('\\Delta=%.2f', mean_val_diff_GB(1)), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fs, 'Color', 'r', 'FontName', 'Times New Roman', 'Interpreter', 'tex');
text(time(2) + 0.2, mean_val_diff_GB(2) + 0.05, sprintf('\\Delta=%.2f', mean_val_diff_GB(2)), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fs, 'Color', 'r', 'FontName', 'Times New Roman', 'Interpreter', 'tex');

% Add p-values and hypothesis test decisions below the bars
% Interest level for Group A
text(time(1) - 0.1, -1, ...
    sprintf('%.1e\n%s', p_int_PC_MC_GA, decision(p_int_PC_MC_GA, alpha_level)), ...
    'HorizontalAlignment', 'center', 'FontSize', fs, 'Color', 'b', 'FontName', 'Times New Roman');
% For Group B
text(time(1) + 0.25, -1, ...
    sprintf('%.1e\n%s', p_int_PC_MC_GB, decision(p_int_PC_MC_GB, alpha_level)), ...
    'HorizontalAlignment', 'center', 'FontSize', fs, 'Color', 'r', 'FontName', 'Times New Roman');

% Interest level for Group A
text(time(2) - 0.25, -1, ...
    sprintf('%.1e\n%s', p_se_PC_MC_GA, decision(p_se_PC_MC_GA, alpha_level)), ...
    'HorizontalAlignment', 'center', 'FontSize', fs, 'Color', 'b', 'FontName', 'Times New Roman');
% For Group B
text(time(2) + 0.25, -1, ...
    sprintf('%.1e\n%s', p_se_PC_MC_GB, decision(p_se_PC_MC_GB, alpha_level)), ...
    'HorizontalAlignment', 'center', 'FontSize', fs, 'Color', 'r', 'FontName', 'Times New Roman');

% Add labels for p-values and hypothesis test decisions
text(0.65, -0.94, ...
    sprintf('P-Val:'), 'HorizontalAlignment', 'center', 'FontSize', fs, 'Color', 'k', 'FontName', 'Times New Roman');
text(0.65, -1.07, ...
    sprintf('Rej H_0:'), 'HorizontalAlignment', 'center', 'FontSize', fs, 'Color', 'k', 'FontName', 'Times New Roman');

hold off