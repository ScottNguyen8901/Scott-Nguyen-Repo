clear
clc
close all

groups = {'A', 'B'};
sheetnames = {'PC', 'MC'};
data_int = cell(2, 2);  % Stores data for Group A & B, PC & MC
data_se = cell(2, 2);   % Stores self-efficacy data for Group A & B, PC & MC
mean_data = cell(2, 2); % Stores mean values for PC & MC

% Loop over groups and sheetnames to read data
for g = 1:2
    for s = 1:2
        data_int{g,s} = readmatrix(['Group ', groups{g}, '\Interest Assessment.xlsx'], 'Sheet', sheetnames{s});
        data_se{g,s} = readmatrix(['Group ', groups{g}, '\Self-Efficacy.xlsx'], 'Sheet', sheetnames{s});
    end
end

% Calculate means and changes
mean_PC_GA = [mean(data_int{1,1}, 2), mean(data_se{1,1}, 2)];
mean_PC_GB = [mean(data_int{2,1}, 2), mean(data_se{2,1}, 2)];
mean_MC_GA = [mean(data_int{1,2}, 2), mean(data_se{1,2}, 2)];
mean_MC_GB = [mean(data_int{2,2}, 2), mean(data_se{2,2}, 2)];

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

% Define parameters for text display
pos = {pos_1, pos_2};
mean_ints = {mean_int_GA, mean_int_GB};
colors = {'b', 'r'};

figure;
subplot(1,2,1)

% Boxplot
h1 = boxplot(mean_int_GA, 'positions', pos_1, 'colors', 'b');
hold on
h2 = boxplot(mean_int_GB, 'positions', pos_2, 'colors', 'r');

% Customize outliers and line widths
set(findobj([h1, h2], 'type', 'line'), 'LineWidth', 2);
set(findobj(h1, 'tag', 'Outliers'), 'MarkerEdgeColor', 'b');
set(gca, 'XTickLabel', x_label, 'XTick', x_tick_values, 'FontSize', 24, 'FontName', 'Times New Roman');

% Display mean values
for j = 1:2
    for i = 1:numel(pos{j})
        mu = mean(mean_ints{j}(:,i));
        text(pos{j}(i), mu + 0.03 * range(ylim), sprintf('\\mu = %.2f', mu), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 11, 'Color', colors{j});
    end
end

% Label and formatting
ylabel('Score', 'FontSize', 24, 'FontName', 'Times New Roman')
axis square;
grid on;

hold off

subplot(1,2,2)

% Boxplot
h3 = boxplot(mean_se_GA, 'positions', pos_1, 'colors', 'b');
hold on
h4 = boxplot(mean_se_GB, 'positions', pos_2, 'colors', 'r');

% Customize outliers and line widths
set(findobj([h3, h4], 'type', 'line'), 'LineWidth', 2);
set(findobj(h3, 'tag', 'Outliers'), 'MarkerEdgeColor', 'b');
set(gca, 'XTickLabel', x_label, 'XTick', x_tick_values, 'FontSize', 24, 'FontName', 'Times New Roman');

% Display mean values
for j = 1:2
    pos_vals = {pos_1, pos_2};
    mean_vals = {mean_se_GA, mean_se_GB};
    colors = {'b', 'r'};
    
    for i = 1:numel(pos_vals{j})
        mu = mean(mean_vals{j}(:,i));
        text(pos_vals{j}(i), mu + 0.02 * range(ylim), sprintf('\\mu = %.2f', mu), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 11, 'Color', colors{j});
    end
end

% Label and formatting
ylabel('Score', 'FontSize', 24, 'FontName', 'Times New Roman')
axis square;
grid on;
legend([h3(1), h4(1)], 'Group A', 'Group B', 'TextColor', 'black', 'Location', 'best');

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

% Add mean value differences on top of the bars for Group A and B
for i = 1:2
    text(time(i) - 0.2, mean_val_diff_GA(i) + 0.05, sprintf('\\Delta=%.2f', mean_val_diff_GA(i)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fs, 'Color', 'b', 'FontName', 'Times New Roman', 'Interpreter', 'tex');
    text(time(i) + 0.2, mean_val_diff_GB(i) + 0.05, sprintf('\\Delta=%.2f', mean_val_diff_GB(i)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', fs, 'Color', 'r', 'FontName', 'Times New Roman', 'Interpreter', 'tex');
end

% Add p-values and hypothesis test decisions below the bars for Group A and B
for i = 1:2
    text(time(i) - 0.1, -1, sprintf('%.1e\n%s', eval(sprintf('p_int_PC_MC_GA')), decision(eval(sprintf('p_int_PC_MC_GA')), alpha_level)), 'HorizontalAlignment', 'center', 'FontSize', fs, 'Color', 'b', 'FontName', 'Times New Roman');
    text(time(i) + 0.25, -1, sprintf('%.1e\n%s', eval(sprintf('p_int_PC_MC_GB')), decision(eval(sprintf('p_int_PC_MC_GB')), alpha_level)), 'HorizontalAlignment', 'center', 'FontSize', fs, 'Color', 'r', 'FontName', 'Times New Roman');
end

% Add labels for p-values and hypothesis test decisions
text(0.65, -0.94, 'P-Val:', 'HorizontalAlignment', 'center', 'FontSize', fs, 'Color', 'k', 'FontName', 'Times New Roman');
text(0.65, -1.07, 'Rej H_0:', 'HorizontalAlignment', 'center', 'FontSize', fs, 'Color', 'k', 'FontName', 'Times New Roman');

hold off