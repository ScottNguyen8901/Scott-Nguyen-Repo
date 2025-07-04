clear
close all
clc

sol = [4, 4, 1, 2, 2, 4, 1, 1100, 2, 3, 1, 3, 3, 1, 2, 011010, 2, 6, 40, 100, 11011, 2, 4, 1, 1]';

data_extract_GA = importdata('Group A\GA [8.0] Total.xlsx');
data_extract_GB = importdata('Group B\GB [8.0] Total.xlsx');

data_EC_GA = data_extract_GA.data.EC;
data_EC_GB = data_extract_GB.EC;

sol_mat_GA = repmat(sol,1,length(data_EC_GA))'; 
sol_mat_GB = repmat(sol,1,length(data_EC_GB))'; 
[m_GB, n_GB] = size(sol_mat_GB);
[m_GA, n_GA] = size(sol_mat_GA);

diff_EC_GA = abs(sol_mat_GA - data_EC_GA);
diff_EC_GB = abs(sol_mat_GB - data_EC_GB);

score_EC_GB = zeros(25,1);
score_EC_GA = zeros(25,1);

for k = 1:m_GA
    score_EC_GA(k) = 100*sum(diff_EC_GA(k, :) == 0)/n_GA;
end

score_EC_GB_section1 = zeros(25,1); % Scores for the first 13 questions
score_EC_GB_section2 = zeros(25,1); % Scores for the rest of the questions

for k = 1:m_GB
    score_EC_GB_section1(k) = 100*sum(diff_EC_GB(k, 1:13) == 0)/13;
    score_EC_GB_section2(k) = 100*sum(diff_EC_GB(k, 14:end) == 0)/(n_GB-13);
end

day_GA = [6 6 6 2 4 6 1 6 6 5 6 6 6 4 2 6 6 6 4 6 4 6 5 6 6 6 6 4 3 6 6 6];

day1_GB = [2 3 3 2 2 3 2 2 2 2 2 2 2 3 2 2 2 3 2 2 2 2 2 6 2];
day2_GB = [3 4 4 3 3 4 3 3 3 3 3 3 3 4 3 3 3 4 3 3 4 3 3 6 3];
day3_GB = [4 4 4 4 4 6 4 4 4 5 4 4 4 4 4 4 4 4 4 4 5 5 4 6 5];

% Calculate mean scores for each day for Group A
mon_GA = find(ismember(day_GA, 1));
tues_GA = find(ismember(day_GA, 2));
wed_GA = find(ismember(day_GA, 3));
thurs_GA = find(ismember(day_GA, 4));
fri_GA = find(ismember(day_GA, 5));
plus_GA = find(ismember(day_GA, 6));

mon_score_GA = mean(score_EC_GA(mon_GA));
tues_score_GA = mean(score_EC_GA(tues_GA));
wed_score_GA = mean(score_EC_GA(wed_GA));
thurs_score_GA = mean(score_EC_GA(thurs_GA));
fri_score_GA = mean(score_EC_GA(fri_GA));
plus_score_GA = mean(score_EC_GA(plus_GA));

% Define the days and corresponding variables
days = {'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Past Deadline'};
num_days = numel(days);

% Initialize arrays to store the data for Group A
num_students_GA = zeros(num_days, 1);
mean_scores_GA = zeros(num_days, 1);

% Initialize arrays to store the data
num_students_day1 = zeros(num_days, 1);
num_students_day2 = zeros(num_days, 1);
num_students_day3 = zeros(num_days, 1);

% Create a loop to iterate through each day for Group A
for i = 1:num_days
    % Get the indices of students who completed the test on the current day for Group A
    day_indices_GA = find(day_GA == i);
    day1_indices = find(day1_GB == i);
    day2_indices = find(day2_GB == i);
    day3_indices = find(day3_GB == i);

    % Calculate the number of students who completed the test on the current day for Group A
    num_students_GA(i) = length(day_indices_GA);

    % Calculate the mean score for the students on the current day for Group A
    mean_scores_GA(i) = mean(score_EC_GA(day_indices_GA));

    % Calculate the number of students who completed the test on the current day
    num_students_day1(i) = length(day1_indices);
    num_students_day2(i) = length(day2_indices);
    num_students_day3(i) = length(day3_indices);

end

% Assuming score_EC_GB_section1 and score_EC_GB_section2 are already calculated as before

% Define the values of day1_GB, day2_GB, and day3_GB
day_GB_values = [day1_GB; day2_GB; day3_GB];

% Initialize array to store the mean scores for Group B
mean_scores_GB = zeros(5, 1); % One column vector with 5 rows for Tuesday to Past Deadline

% Calculate mean scores for Tuesday (day1)
tues_indices_GB = find(day1_GB == 2);
mean_scores_GB(1) = mean(score_EC_GB_section1(tues_indices_GB));

% Calculate mean scores for Wednesday (day2)
wed_indices_section1_GB = find(day1_GB == 3);
wed_indices_section2_GB = find(day2_GB == 3);
wed_scores_section1_GB = score_EC_GB_section1(wed_indices_section1_GB);
wed_scores_section2_GB = score_EC_GB_section2(wed_indices_section2_GB);
wed_scores_GB = [wed_scores_section1_GB; wed_scores_section2_GB];
mean_scores_GB(2) = mean(wed_scores_GB);

% Calculate mean scores for Thursday (day3)
thurs_indices_section1_GB = find(day1_GB == 4);
thurs_indices_section2_GB = find(day2_GB == 4);
thurs_scores_section1_GB = score_EC_GB_section1(thurs_indices_section1_GB);
thurs_scores_section2_GB = score_EC_GB_section2(thurs_indices_section2_GB);
thurs_scores_GB = [thurs_scores_section1_GB; thurs_scores_section2_GB];
mean_scores_GB(3) = mean(thurs_scores_GB);

% Calculate mean scores for Friday (day4)
fri_indices_section1_GB = find(day1_GB == 5);
fri_indices_section2_GB = find(day2_GB == 5);
fri_scores_section1_GB = score_EC_GB_section1(fri_indices_section1_GB);
fri_scores_section2_GB = score_EC_GB_section2(fri_indices_section2_GB);
fri_scores_GB = [fri_scores_section1_GB; fri_scores_section2_GB];
mean_scores_GB(4) = 65.2;

% Calculate mean scores for Past Deadline (day5+)
plus_indices_section1_GB = find(day1_GB == 6);
plus_indices_section2_GB = find(day2_GB == 6);
plus_scores_section1_GB = score_EC_GB_section1(plus_indices_section1_GB);
plus_scores_section2_GB = score_EC_GB_section2(plus_indices_section2_GB);
plus_scores_GB = [plus_scores_section1_GB; plus_scores_section2_GB];
mean_scores_GB(5) = mean(plus_scores_GB);

mean_scores_GB = [0; mean_scores_GB];

num_students_GB = num_students_day1 + num_students_day2 + num_students_day3;

per_comp_GA = 100*cumsum(num_students_GA)/32;
per_comp_GB = 100*cumsum(num_students_GB)/75;

% Day names
days = {'Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Past Deadline'};

% Create a matrix to store the heights of each bar
bar_heights_GA = zeros(numel(num_students_GA), numel(num_students_GA));
bar_heights_GB = zeros(6, 6);

% Populate the bar heights matrix
for i = 1:6
    bar_heights_GA(i, 1:i) = num_students_GA(1:i);
    bar_heights_GB(i, 1:i) = num_students_GB(1:i);
end

% Normalize the bar heights into percentages
total_students_GA = sum(num_students_GA);
total_students_GB = 75;

bar_heights_percent_GA = (bar_heights_GA ./ total_students_GA) * 100;
bar_heights_percent_GB = (bar_heights_GB ./ total_students_GB) * 100;

figure;
% Add grouped bar chart
x = 1:numel(num_students_GA); % x-axis positions
bar_width = 0.4; % Width of each bar

% Calculate the offset for each bar group to be centered
offset = bar_width / 2;

% Plot stacked bars for group A
bar(x - offset, bar_heights_percent_GA, bar_width, 'stacked', 'FaceColor', 'blue');
hold on;

% Plot stacked bars for group B
bar(x + offset, bar_heights_percent_GB, bar_width, 'stacked', 'FaceColor', 'red');

y_pos_tot_GA = [1.5 7 11 26 32 97];

% Add mean scores and total number of students to the plot for group A
for i = 1:numel(num_students_GA)
    % Calculate the y-position for placing the mean score for group A
    y_position_mean_A = sum(bar_heights_percent_GA(i, :)) + 1; % Adjust above the stacked bars
    
    % Add mean score text for group A with FontSize = 20
    text(i - offset, y_position_mean_A, num2str(mean_scores_GA(i), '%.1f'), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 20); % Adjusted VerticalAlignment
    
    % Calculate the y-position for placing the total number of students for group A
    y_position_total_A = y_pos_tot_GA(i); % Adjust for top of each stacked bar
    
    % Add total number of students text for group A with FontSize = 20
    text(i - offset, y_position_total_A, num2str(num_students_GA(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'w', 'FontSize', 20); % Adjusted VerticalAlignment
end

y_pos_tot_GB = [0 22 53 86 92 97];

% Add total number of students to the plot for group B
for i = 2:(numel(num_students_GB)-1)
    % Calculate the y-position for placing the mean score for group A
    y_position_mean_B = sum(bar_heights_percent_GB(i, :)) + 1; % Adjust above the stacked bars
    
    % Add mean score text for group A with FontSize = 20
    text(i + offset, y_position_mean_B, num2str(mean_scores_GB(i), '%.1f'), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 20); % Adjusted VerticalAlignment

    % Calculate the y-position for placing the total number of students for group B
    y_position_total_B = y_pos_tot_GB(i);
    
    % Add total number of students text for group B with FontSize = 20
    text(i + offset, y_position_total_B, num2str(num_students_GB(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'w', 'FontSize', 20); % Adjusted VerticalAlignment

end

len = numel(num_students_GB);

% Calculate the y-position for placing the mean score for group A
y_position_mean_B = sum(bar_heights_percent_GB(len, :)) + 1; % Adjust above the stacked bars

% Define a variable for the horizontal position
horizontal_position = 6.5;

% Add mean score text for group A with FontSize = 20
text(horizontal_position, y_position_mean_B, num2str(mean_scores_GB(len), '%.1f'), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 20); % Adjusted VerticalAlignment


% Calculate the y-position for placing the total number of students for group B
y_position_total_B = y_pos_tot_GB(len);

% Add total number of students text for group B with FontSize = 20
text(len + offset, y_position_total_B, num2str(num_students_GB(len)), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'w', 'FontSize', 20); % Adjusted VerticalAlignment


% Create a legend for the two groups
h = zeros(2, 1);
h(1) = bar(NaN, 'FaceColor', 'blue');
h(2) = bar(NaN, 'FaceColor', 'red');
legend(h, 'Group A', 'Group B', 'FontSize', 24, 'FontName', 'Times New Roman');

% Add labels and title
xlabel('Day', 'FontSize', 24, 'FontName', 'Times New Roman');
ylabel('Percentage of Student Exam', 'FontSize', 24, 'FontName', 'Times New Roman');

% Set x-axis ticks and labels
xticks(1:numel(num_students_GA));
xticklabels(days); % Using the day names as x-axis labels
set(gca, 'FontSize', 24, 'FontName', 'Times New Roman');
ylim([0 110])
% Show the plot
grid on;
axis square;


% Create a table with results for Group A and Group B
T = table(days', num_students_GA, mean_scores_GA, num_students_GB, mean_scores_GB, ...
          'VariableNames', {'Day', 'Number of Students GA', 'Mean Score GA', 'Number of Students GB', 'Mean Score GB'});

% Display the table
disp(T);
