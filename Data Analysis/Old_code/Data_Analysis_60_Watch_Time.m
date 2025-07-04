clear, clc, close all

data_extract_GA = importdata('Watch Times\GA AE 298 Video Information.xlsx');
data_extract_GB = importdata('Watch Times\GB AE 298 Video Information.xlsx');

data_GA = data_extract_GA.data.VideoAnalytics;
data_GA(:,3) = [];
data_GA(:,5) = [];
data_GA(:,end) = [];

data_GB = data_extract_GB.data.VideoAnalytics;
data_GB(:,3) = [];
data_GB(:,5) = [];
data_GB(:,end) = [];

data_GA_buffer = data_extract_GA.data.VideoAnalyticsBuffer;
data_GA_buffer(:,3) = [];
data_GA_buffer(:,5) = [];
data_GA_buffer(:,end) = [];

data_GB_buffer = data_extract_GB.data.VideoAnalyticsBuffer;
data_GB_buffer(:,3) = [];
data_GB_buffer(:,5) = [];
data_GB_buffer(:,end) = [];

data_GA_final = data_extract_GA.data.VideoAnalyticsFinal;
data_GA_final(:,3) = [];
data_GA_final(:,5) = [];
data_GA_final(:,end) = [];

data_GB_final = data_extract_GB.data.VideoAnalyticsFinal;
data_GB_final(:,3) = [];
data_GB_final(:,5) = [];
data_GB_final(:,end) = [];

mod_1_data_GA = data_GA(1:11,:);
mod_1_data_GA(4,:) = [];
mod_2_data_GA = data_GA(13:21,:);
mod_3_data_GA = data_GA(23:26,:);
mod_4_data_GA = data_GA(28:end,:);

mod_1_data_GB = data_GB(1:11,:);
mod_1_data_GB(4,:) = [];
mod_2_data_GB = data_GB(13:21,:);
mod_3_data_GB = data_GB(23:26,:);
mod_4_data_GB = data_GB(28:end,:);

mod_1_data_GA_buffer = data_GA_buffer(1:11,:);
mod_1_data_GA_buffer(4,:) = [];
mod_2_data_GA_buffer = data_GA_buffer(13:21,:);
mod_3_data_GA_buffer = data_GA_buffer(23:26,:);
mod_4_data_GA_buffer = data_GA_buffer(28:end,:);

mod_1_data_GB_buffer = data_GB_buffer(1:11,:);
mod_1_data_GB_buffer(4,:) = [];
mod_2_data_GB_buffer = data_GB_buffer(13:21,:);
mod_3_data_GB_buffer = data_GB_buffer(23:26,:);
mod_4_data_GB_buffer = data_GB_buffer(28:end,:);

mod_1_data_GB_Final = data_GB_final(1:11,:);
mod_1_data_GB_Final(4,:) = [];
mod_2_data_GB_Final = data_GB_final(13:21,:);
mod_3_data_GB_Final = data_GB_final(23:26,:);
mod_4_data_GB_Final = data_GB_final(28:end,:);

mod_1_data_GA_Final = data_GA_final(1:11,:);
mod_1_data_GA_Final(4,:) = [];
mod_2_data_GA_Final = data_GA_final(13:21,:);
mod_3_data_GA_Final = data_GA_final(23:26,:);
mod_4_data_GA_Final = data_GA_final(28:end,:);

mod_1_video_names = {
    'Introduction to Single Stage Rocket'
    'Why we go to Space'
    'Phases of Rocket Flight'
    'Introduction to Rocket Hardware'
    'Rocket Bodies'
    'Rocket Engines'
    'Recovery Systems'
    'Launch Controller'
    'Avionics'
    'Payload'
};

mod_2_video_names = {
    'Introduction to Rocket Design'
    'Center of Gravity'
    'Center of Pressure'
    'Equilibrium'
    'Low Velocity Stability'
    'High Velocity Stability'
    'Thrust, Weight, and Impulse'
    'Thrust to Weight Ratio'
    'Motor Selection'
};

mod_3_video_names = {
    'Introduction to Rocket Mechanics'
    'Derive and Describe Rocket EOMs'
    'Solving Approximate EOMs for Altitude (No Drag)'
    'Plotting Altitude (Google Sheets)'
};

mod_4_video_names = {
    'Introduction to Analysis'
    'Comparing Different Models Part 1'
    'Comparing Different Models Part 2'
    'Compare Flight Data to Predictions'
    'Discussion with an Astra Engineer'
    'A Discussion with an Aerospace Engineering Student'
    'The Future of Space and Rocketry'
    'Conclusion'
};

video_names = [mod_1_video_names; mod_2_video_names; mod_3_video_names; mod_4_video_names];

len_1 = length(mod_1_data_GB(:,2));
len_2 = length(mod_2_data_GB(:,2));
len_3 = length(mod_3_data_GB(:,2));
len_4 = length(mod_4_data_GB(:,2));

% figure;
% subplot(1,2,1)
% plot(1:len_1, mod_1_data_GA(:,2), 'b-o');
% hold on
% plot(len_1 + 1:len_1 + len_2, mod_2_data_GA(:,2), 'r-o');
% plot(len_1 + len_2 + 1:len_1 + len_2 + len_3, mod_3_data_GA(:,2), 'm-o');
% plot(len_1 + len_2 + len_3 + 1:len_1 + len_2 + len_3 + len_4, mod_4_data_GA(:,2), 'k-o');
% plot(1:len_1, mod_1_data_GB(:,2), 'b--*');
% plot(len_1 + 1:len_1 + len_2, mod_2_data_GB(:,2), 'r--*');
% plot(len_1 + len_2 + 1:len_1 + len_2 + len_3, mod_3_data_GB(:,2), 'm--*');
% plot(len_1 + len_2 + len_3 + 1:len_1 + len_2 + len_3 + len_4, mod_4_data_GB(:,2), 'k--*');
% grid on
% 
% % Set x-axis labels
% xticks(1:numel(video_names));
% xticklabels(video_names);
% xtickangle(45); % Rotate x-axis labels for better readability
% 
% % Set axis labels and title
% xlabel('Video Names');
% ylabel('Views/Student');
% title('Views/Student for Videos');
% 
% subplot(1,2,2)
% plot(1:len_1, mod_1_data_GA(:,4), 'b-o');
% hold on
% plot(len_1 + 1:len_1 + len_2, mod_2_data_GA(:,4), 'r-o');
% plot(len_1 + len_2 + 1:len_1 + len_2 + len_3, mod_3_data_GA(:,4), 'm-o');
% plot(len_1 + len_2 + len_3 + 1:len_1 + len_2 + len_3 + len_4, mod_4_data_GA(:,4), 'k-o');
% plot(1:len_1, mod_1_data_GB(:,4), 'b--*');
% plot(len_1 + 1:len_1 + len_2, mod_2_data_GB(:,4), 'r--*');
% plot(len_1 + len_2 + 1:len_1 + len_2 + len_3, mod_3_data_GB(:,4), 'm--*');
% plot(len_1 + len_2 + len_3 + 1:len_1 + len_2 + len_3 + len_4, mod_4_data_GB(:,4), 'k--*');
% grid on
% legend( ...
%     'Group A Introduction and Rocket Hardware', ...
%     'Group A Fundamentals of Rocketry', ...
%     'Group A Modeling Rocket Dynamics', ...
%     'Group A Analysis', ...
%     'Group B Introduction and Rocket Hardware', ...
%     'Group B Fundamentals of Rocketry', ... 
%     'Group B Modeling Rocket Dynamics', ...
%     'Group B Analysis')
% 
% % Set x-axis labels
% xticks(1:numel(video_names));
% xticklabels(video_names);
% xtickangle(45); % Rotate x-axis labels for better readability
% 
% % Set axis labels and title
% xlabel('Video Names');
% ylabel('View Duration (%)');
% title('Percent View Duration for Videos');
% 
% sgtitle('Video Comparison No Buffer');
% 
% figure(2)
% subplot(1,2,1)
% plot(1:len_1, mod_1_data_GA_buffer(:,2), 'b-o');
% hold on
% plot(len_1 + 1:len_1 + len_2, mod_2_data_GA_buffer(:,2), 'r-o');
% plot(len_1 + len_2 + 1:len_1 + len_2 + len_3, mod_3_data_GA_buffer(:,2), 'm-o');
% plot(len_1 + len_2 + len_3 + 1:len_1 + len_2 + len_3 + len_4, mod_4_data_GA_buffer(:,2), 'k-o');
% plot(1:len_1, mod_1_data_GB(:,2), 'b--*');
% plot(len_1 + 1:len_1 + len_2, mod_2_data_GB_buffer(:,2), 'r--*');
% plot(len_1 + len_2 + 1:len_1 + len_2 + len_3, mod_3_data_GB_buffer(:,2), 'm--*');
% plot(len_1 + len_2 + len_3 + 1:len_1 + len_2 + len_3 + len_4, mod_4_data_GB_buffer(:,2), 'k--*');
% grid on
% 
% % Set x-axis labels
% xticks(1:numel(video_names));
% xticklabels(video_names);
% xtickangle(45); % Rotate x-axis labels for better readability
% 
% % Set axis labels and title
% xlabel('Video Names');
% ylabel('Views/Student');
% title('Views/Student for Videos');
% 
% subplot(1,2,2)
% plot(1:len_1, mod_1_data_GA(:,4), 'b-o');
% hold on
% plot(len_1 + 1:len_1 + len_2, mod_2_data_GA_buffer(:,4), 'r-o');
% plot(len_1 + len_2 + 1:len_1 + len_2 + len_3, mod_3_data_GA_buffer(:,4), 'm-o');
% plot(len_1 + len_2 + len_3 + 1:len_1 + len_2 + len_3 + len_4, mod_4_data_GA_buffer(:,4), 'k-o');
% plot(1:len_1, mod_1_data_GB_buffer(:,4), 'b--*');
% plot(len_1 + 1:len_1 + len_2, mod_2_data_GB_buffer(:,4), 'r--*');
% plot(len_1 + len_2 + 1:len_1 + len_2 + len_3, mod_3_data_GB_buffer(:,4), 'm--*');
% plot(len_1 + len_2 + len_3 + 1:len_1 + len_2 + len_3 + len_4, mod_4_data_GB_buffer(:,4), 'k--*');
% grid on
% legend( ...
%     'Group A Introduction and Rocket Hardware', ...
%     'Group A Fundamentals of Rocketry', ...
%     'Group A Modeling Rocket Dynamics', ...
%     'Group A Analysis', ...
%     'Group B Introduction and Rocket Hardware', ...
%     'Group B Fundamentals of Rocketry', ... 
%     'Group B Modeling Rocket Dynamics', ...
%     'Group B Analysis')
% 
% % Set x-axis labels
% xticks(1:numel(video_names));
% xticklabels(video_names);
% xtickangle(45); % Rotate x-axis labels for better readability
% 
% % Set axis labels and title
% xlabel('Video Names');
% ylabel('View Duration (%)');
% title('Percent View Duration for Videos');
% 
% sgtitle('Video Comparison with Buffer');
% 
% figure;
% subplot(1,2,1)
% plot(1:len_1, 32*mod_1_data_GA(:,2), 'b-o');
% hold on
% plot(len_1 + 1:len_1 + len_2, 32*mod_2_data_GA(:,2), 'r-o');
% plot(len_1 + len_2 + 1:len_1 + len_2 + len_3, 32*mod_3_data_GA(:,2), 'm-o');
% plot(len_1 + len_2 + len_3 + 1:len_1 + len_2 + len_3 + len_4, 32*mod_4_data_GA(:,2), 'k-o');
% plot(1:len_1, 25*mod_1_data_GB(:,2), 'b--*');
% plot(len_1 + 1:len_1 + len_2, 25*mod_2_data_GB(:,2), 'r--*');
% plot(len_1 + len_2 + 1:len_1 + len_2 + len_3, 25*mod_3_data_GB(:,2), 'm--*');
% plot(len_1 + len_2 + len_3 + 1:len_1 + len_2 + len_3 + len_4, 25*mod_4_data_GB(:,2), 'k--*');
% grid on
% 
% % Set x-axis labels
% xticks(1:numel(video_names));
% xticklabels(video_names);
% xtickangle(45); % Rotate x-axis labels for better readability
% 
% % Set axis labels and title
% xlabel('Video Names');
% ylabel('Views');
% title('Views for Videos');
% 
% subplot(1,2,2)
% plot(1:len_1, mod_1_data_GA(:,3), 'b-o');
% hold on
% plot(len_1 + 1:len_1 + len_2, mod_2_data_GA(:,3), 'r-o');
% plot(len_1 + len_2 + 1:len_1 + len_2 + len_3, mod_3_data_GA(:,3), 'm-o');
% plot(len_1 + len_2 + len_3 + 1:len_1 + len_2 + len_3 + len_4, mod_4_data_GA(:,3), 'k-o'); 
% plot(1:len_1, mod_1_data_GB(:,4), 'b--*');
% plot(len_1 + 1:len_1 + len_2, mod_2_data_GB(:,3), 'r--*');
% plot(len_1 + len_2 + 1:len_1 + len_2 + len_3, mod_3_data_GB(:,3), 'm--*');
% plot(len_1 + len_2 + len_3 + 1:len_1 + len_2 + len_3 + len_4, mod_4_data_GB(:,3), 'k--*');
% grid on
% legend( ...
%     'Group A Introduction and Rocket Hardware', ...
%     'Group A Fundamentals of Rocketry', ...
%     'Group A Modeling Rocket Dynamics', ...
%     'Group A Analysis', ...
%     'Group B Introduction and Rocket Hardware', ...
%     'Group B Fundamentals of Rocketry', ... 
%     'Group B Modeling Rocket Dynamics', ...
%     'Group B Analysis')
% 
% % Set x-axis labels
% xticks(1:numel(video_names));
% xticklabels(video_names);
% xtickangle(45); % Rotate x-axis labels for better readability
% 
% % Set axis labels and title
% xlabel('Video Names');
% ylabel('View Duration (s)');
% title('View Duration for Videos');
% 
% sgtitle('Video Comparison No Buffer');

% Assuming len_1, len_2, len_3, len_4, mod_1_data_GA_buffer, etc. are defined elsewhere

% Define the number of videos (assuming this is known)
num_videos = len_1 + len_2 + len_3 + len_4;
module_boundaries = [1, len_1, len_1 + len_2, len_1 + len_2 + len_3, num_videos];

% Calculate the midpoints for x-axis labels
midpoints = (module_boundaries(1:end-1) + module_boundaries(2:end)) / 2;

% Start a new figure
figure;
h1 = plot(1:len_1, 32*mod_1_data_GA_buffer(:,2), 'bo-', 'LineWidth', 2);
hold on;
plot(len_1 + 1:len_1 + len_2, 32*mod_2_data_GA_buffer(:,2), 'bo-', 'LineWidth', 2);
plot(len_1 + len_2 + 1:len_1 + len_2 + len_3, 32*mod_3_data_GA_buffer(:,2), 'bo-', 'LineWidth', 2);
plot(len_1 + len_2 + len_3 + 1:num_videos, 32*mod_4_data_GA_buffer(:,2), 'bo-', 'LineWidth', 2);
h2 = plot(1:len_1, mod_1_data_GB_buffer(:,1), 'r--*', 'LineWidth', 2);
plot(len_1 + 1:len_1 + len_2, 25*mod_2_data_GB_buffer(:,2), 'r--*', 'LineWidth', 2);
plot(len_1 + len_2 + 1:len_1 + len_2 + len_3, 25*mod_3_data_GB_buffer(:,2), 'r--*', 'LineWidth', 2);
plot(len_1 + len_2 + len_3 + 1:num_videos, 25*mod_4_data_GB_buffer(:,2), 'r--*', 'LineWidth', 2);
grid on;
axis square;

% Set x-axis labels in increments of 5
xticks(0:5:num_videos);
xlim([0 num_videos]); % Set the limits of the x-axis

% Customize x and y axis ticks
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);

% Add a single horizontal bar at y=65
line(xlim, [65, 65], 'Color', 'k', 'LineWidth', 2);

% Add small vertical bars at the boundaries of each module
for i = 2:length(module_boundaries)-1
    line([module_boundaries(i), module_boundaries(i)], [65 - 2, 65 + 2], 'Color', 'k', 'LineWidth', 2); % 10 units above and below the horizontal line
end

% Print module labels above the horizontal bar
for i = 1:length(midpoints)
    text(midpoints(i), 65 + 2, sprintf('M%d', i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontName', 'Times New Roman', 'FontSize', 24); % Adjust the y-value as needed
end

% Create a legend for the first plot of each group
legend([h1, h2], {'Group A', 'Group B'}, 'FontName', 'Times New Roman', 'FontSize', 24);

% Set y-axis label
xlabel('Video Number', 'FontName', 'Times New Roman', 'FontSize', 24);
ylabel('Views', 'FontName', 'Times New Roman', 'FontSize', 24);

figure
% Subplot 1
subplot(1, 2, 1);
% Plot Group A data and set the legend 'IconDisplayStyle' to 'off' for lines
h1 = plot(1:len_1, mod_1_data_GA(:,3), 'bo-', 'LineWidth', 2); % Group A solid black with open circle
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','on'); % Include in legend
hold on;
% For the rest of Group A data, exclude lines from legend
h = plot(len_1 + 1:len_1 + len_2, mod_2_data_GA_buffer(:,3), 'bo-', 'LineWidth', 2);
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend
h = plot(len_1 + len_2 + 1:len_1 + len_2 + len_3, mod_3_data_GA_buffer(:,3), 'bo-', 'LineWidth', 2);
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend
h = plot(len_1 + len_2 + len_3 + 1:num_videos, mod_4_data_GA_buffer(:,3), 'bo-', 'LineWidth', 2);
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend

% Plot Group B data and set the legend 'IconDisplayStyle' to 'off' for lines
h2 = plot(1:len_1, mod_1_data_GB_buffer(:,3), 'r--*', 'LineWidth', 2); % Group B dashed black with filled star
set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','on'); % Include in legend
% For the rest of Group B data, exclude lines from legend
h = plot(len_1 + 1:len_1 + len_2, mod_2_data_GB_buffer(:,3), 'r--*', 'LineWidth', 2);
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend
h = plot(len_1 + len_2 + 1:len_1 + len_2 + len_3, mod_3_data_GB_buffer(:,3), 'r--*', 'LineWidth', 2);
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend
h = plot(len_1 + len_2 + len_3 + 1:num_videos, mod_4_data_GB_buffer(:,3), 'r--*', 'LineWidth', 2);
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend

grid on;
axis square;

% Set x-axis labels in increments of 5
xticks(0:5:num_videos);
xlim([0 num_videos]); % Set the limits of the x-axis

% Customize x and y axis ticks
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);

% Add a single horizontal bar at y=375 and exclude it from the legend
hLine = line(xlim, [375, 375], 'Color', 'k', 'LineWidth', 2);
set(get(get(hLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% Add small vertical bars at the boundaries of each module and exclude them from the legend
for i = 2:length(module_boundaries)-1
    hLine = line([module_boundaries(i), module_boundaries(i)], [375 - 10, 375 + 10], 'Color', 'k', 'LineWidth', 2); % 10 units above and below the horizontal line
    set(get(get(hLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end

% Print module labels above the horizontal bar
for i = 1:length(midpoints)
    text(midpoints(i), 375 + 5, sprintf('M%d', i), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 24, 'FontName', 'Times New Roman'); % Adjust the y-value as needed
end

% Set y-axis label
xlabel('Video Number', 'FontName', 'Times New Roman', 'FontSize', 24);
ylabel('View Duration (s)', 'FontName', 'Times New Roman', 'FontSize', 24);

% Subplot 2
subplot(1, 2, 2);
% Plot Group A data and set the legend 'IconDisplayStyle' to 'off' for lines after the first
h1 = plot(1:len_1, mod_1_data_GA(:,4), 'bo-', 'LineWidth', 2); % Group A solid blue with open circle
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','on'); % Include in legend
hold on;
h = plot(len_1 + 1:len_1 + len_2, mod_2_data_GA(:,4), 'bo-', 'LineWidth', 2);
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend
h = plot(len_1 + len_2 + 1:len_1 + len_2 + len_3, mod_3_data_GA(:,4), 'bo-', 'LineWidth', 2);
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend
h = plot(len_1 + len_2 + len_3 + 1:len_1 + len_2 + len_3 + len_4, mod_4_data_GA(:,4), 'bo-', 'LineWidth', 2);
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend

% Plot Group B data and set the legend 'IconDisplayStyle' to 'off' for lines after the first
h2 = plot(1:len_1, mod_1_data_GB(:,4), 'r--*', 'LineWidth', 2); % Group B dashed blue with filled star
set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','on'); % Include in legend
h = plot(len_1 + 1:len_1 + len_2, mod_2_data_GB(:,4), 'r--*', 'LineWidth', 2);
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend
h = plot(len_1 + len_2 + 1:len_1 + len_2 + len_3, mod_3_data_GB(:,4), 'r--*', 'LineWidth', 2);
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend
h = plot(len_1 + len_2 + len_3 + 1:len_1 + len_2 + len_3 + len_4, mod_4_data_GB(:,4), 'r--*', 'LineWidth', 2);
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend

grid on;
axis square;

% Set x-axis labels in increments of 5
xticks(0:5:len_1 + len_2 + len_3 + len_4);
xlim([0 len_1 + len_2 + len_3 + len_4]); % Set the limits of the x-axis

% Customize x and y axis ticks
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);

% Add a single horizontal bar at y=375 and exclude it from the legend
hLine = line(xlim, [100, 100], 'Color', 'k', 'LineWidth', 2);
set(get(get(hLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% Add small vertical bars at the boundaries of each module and exclude them from the legend
for i = 2:length(module_boundaries)-1
    hLine = line([module_boundaries(i), module_boundaries(i)], [100 - 5, 100 + 5], 'Color', 'k', 'LineWidth', 2); % Span the y-axis
    set(get(get(hLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end

% Print module labels above the horizontal bar
for i = 1:length(midpoints)
    text(midpoints(i), 100 + 2, sprintf('M%d', i), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 24, 'FontName', 'Times New Roman'); % Adjust the y-value as needed
end

% Create a legend for the first plot of each group
legend([h1, h2], {'Group A', 'Group B'}, 'FontName', 'Times New Roman', 'FontSize', 24);

% Set y-axis label
xlabel('Video Number', 'FontName', 'Times New Roman', 'FontSize', 24);
ylabel('View Duration (%)', 'FontName', 'Times New Roman', 'FontSize', 24);
ylim([0 110])
% figure;
% subplot(1,2,1)
% plot(1:len_1, mod_1_data_GA_Final(:,2), 'b:o');
% hold on
% plot(len_1 + 1:len_1 + len_2, mod_2_data_GA_Final(:,2), 'r:o');
% plot(len_1 + len_2 + 1:len_1 + len_2 + len_3, mod_3_data_GA_Final(:,2), 'm:o');
% plot(len_1 + len_2 + len_3 + 1:len_1 + len_2 + len_3 + len_4, mod_4_data_GA_Final(:,2), 'k:o');
% plot(1:len_1, mod_1_data_GB_Final(:,2), 'b-.*');
% plot(len_1 + 1:len_1 + len_2, mod_2_data_GB_Final(:,2), 'r-.*');
% plot(len_1 + len_2 + 1:len_1 + len_2 + len_3, mod_3_data_GB_Final(:,2), 'm-.*');
% plot(len_1 + len_2 + len_3 + 1:len_1 + len_2 + len_3 + len_4, mod_4_data_GB_Final(:,2), 'k-.*');
% grid on
% 
% % Set x-axis labels
% xticks(1:numel(video_names));
% xticklabels(video_names);
% xtickangle(45); % Rotate x-axis labels for better readability
% 
% % Set axis labels and title
% xlabel('Video Names', 'FontName', 'Times New Roman', 'FontSize', 12);
% ylabel('Views/Student', 'FontName', 'Times New Roman', 'FontSize', 12);
% title('Views/Student for Videos', 'FontName', 'Times New Roman', 'FontSize', 12);
% 
% subplot(1,2,2)
% plot(1:len_1, mod_1_data_GA_Final(:,4), 'b:o');
% hold on
% plot(len_1 + 1:len_1 + len_2, mod_2_data_GA_Final(:,4), 'r:o');
% plot(len_1 + len_2 + 1:len_1 + len_2 + len_3, mod_3_data_GA_Final(:,4), 'm:o');
% plot(len_1 + len_2 + len_3 + 1:len_1 + len_2 + len_3 + len_4, mod_4_data_GA_Final(:,4), 'k:o');
% plot(1:len_1, mod_1_data_GB_Final(:,4), 'b-.*');
% plot(len_1 + 1:len_1 + len_2, mod_2_data_GB_Final(:,4), 'r-.*');
% plot(len_1 + len_2 + 1:len_1 + len_2 + len_3, mod_3_data_GB_Final(:,4), 'm-.*');
% plot(len_1 + len_2 + len_3 + 1:len_1 + len_2 + len_3 + len_4, mod_4_data_GB_Final(:,4), 'k-.*');
% grid on
% legend( ...
%     'Group A Introduction and Rocket Hardware', ...
%     'Group A Fundamentals of Rocketry', ...
%     'Group A Modeling Rocket Dynamics', ...
%     'Group A Analysis', ...
%     'Group B Introduction and Rocket Hardware', ...
%     'Group B Fundamentals of Rocketry', ... 
%     'Group B Modeling Rocket Dynamics', ...
%     'Group B Analysis', ...
%     'FontName', 'Times New Roman', 'FontSize', 12)
% 
% % Set x-axis labels
% xticks(1:numel(video_names));
% xticklabels(video_names);
% xtickangle(45); % Rotate x-axis labels for better readability
% 
% % Set axis labels and title
% xlabel('Video Names', 'FontName', 'Times New Roman', 'FontSize', 12);
% ylabel('View Duration (%)', 'FontName', 'Times New Roman', 'FontSize', 12);
% title('Views/Student for Videos', 'FontName', 'Times New Roman', 'FontSize', 12);
% 
% sgtitle('Finals Week Video Comparison', 'FontName', 'Times New Roman', 'FontSize', 12);

% Start a new figure
figure;

% Subplot 1
subplot(1, 2, 1);
plot(1:len_1, 32*mod_1_data_GA_Final(:,2), 'bo-', 'LineWidth', 2);
hold on;
plot(len_1 + 1:len_1 + len_2, 32*mod_2_data_GA_Final(:,2), 'bo-', 'LineWidth', 2);
plot(len_1 + len_2 + 1:len_1 + len_2 + len_3, 32*mod_3_data_GA_Final(:,2), 'bo-', 'LineWidth', 2);
plot(len_1 + len_2 + len_3 + 1:num_videos, 32*mod_4_data_GA_Final(:,2), 'bo-', 'LineWidth', 2);
plot(1:len_1, 25*mod_1_data_GB_Final(:,2), 'r--*', 'LineWidth', 2);
plot(len_1 + 1:len_1 + len_2, 25*mod_2_data_GB_Final(:,2), 'r--*', 'LineWidth', 2);
plot(len_1 + len_2 + 1:len_1 + len_2 + len_3, 25*mod_3_data_GB_Final(:,2), 'r--*', 'LineWidth', 2);
plot(len_1 + len_2 + len_3 + 1:num_videos, 25*mod_4_data_GB_Final(:,2), 'r--*', 'LineWidth', 2);
grid on;
axis square;

% Set x-axis labels in increments of 5
xticks(0:5:num_videos);
xlim([0 num_videos]); % Set the limits of the x-axis

% Customize x and y axis ticks
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);

% Add a single horizontal bar at y=17
line(xlim, [17, 17], 'Color', 'k', 'LineWidth', 2);

% Add small vertical bars at the boundaries of each module with reduced height
for i = 2:length(module_boundaries)-1
    line([module_boundaries(i), module_boundaries(i)], [17 - 0.5, 17 + 0.5], 'Color', 'k', 'LineWidth', 2); % 1 unit above and below the horizontal line
end

% Print module labels above the horizontal bar
for i = 1:length(midpoints)
    text(midpoints(i), 17, sprintf('M%d', i), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 24, 'FontName', 'Times New Roman'); % Adjust the y-value as needed
end

% Set y-axis label
xlabel('Video Number', 'FontName', 'Times New Roman', 'FontSize', 24);
ylabel('Views', 'FontName', 'Times New Roman', 'FontSize', 24);

% Subplot 2
subplot(1, 2, 2);
h1 = plot(1:len_1, mod_1_data_GA_Final(:,3), 'bo-', 'LineWidth', 2); % Group A solid black with open circle
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','on'); % Include in legend
hold on;
% For the rest of Group A data, exclude lines from legend
h = plot(len_1 + 1:len_1 + len_2, mod_2_data_GA_Final(:,3), 'bo-', 'LineWidth', 2);
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend
h = plot(len_1 + len_2 + 1:len_1 + len_2 + len_3, mod_3_data_GA_Final(:,3), 'bo-', 'LineWidth', 2);
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend
h = plot(len_1 + len_2 + len_3 + 1:num_videos, mod_4_data_GA_Final(:,3), 'bo-', 'LineWidth', 2);
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend

% Plot Group B data and set the legend 'IconDisplayStyle' to 'off' for lines
h2 = plot(1:len_1, mod_1_data_GB_Final(:,3), 'r--*', 'LineWidth', 2); % Group B dashed black with filled star
set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','on'); % Include in legend
% For the rest of Group B data, exclude lines from legend
h = plot(len_1 + 1:len_1 + len_2, mod_2_data_GB_Final(:,3), 'r--*', 'LineWidth', 2);
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend
h = plot(len_1 + len_2 + 1:len_1 + len_2 + len_3, mod_3_data_GB_Final(:,3), 'r--*', 'LineWidth', 2);
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend
h = plot(len_1 + len_2 + len_3 + 1:num_videos, mod_4_data_GB_Final(:,3), 'r--*', 'LineWidth', 2);
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude from legend

grid on;
axis square;

% Set x-axis labels in increments of 5
xticks(0:5:num_videos);
xlim([0 num_videos]); % Set the limits of the x-axis

% Customize x and y axis ticks
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);

% Add a single horizontal bar at y=375 and exclude it from the legend
hLine = line(xlim, [320, 320], 'Color', 'k', 'LineWidth', 2);
set(get(get(hLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% Add small vertical bars at the boundaries of each module and exclude them from the legend
for i = 2:length(module_boundaries)-1
    hLine = line([module_boundaries(i), module_boundaries(i)], [320 - 10, 320 + 10], 'Color', 'k', 'LineWidth', 2); % 10 units above and below the horizontal line
    set(get(get(hLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end

% Print module labels above the horizontal bar
for i = 1:length(midpoints)
    text(midpoints(i), 320, sprintf('M%d', i), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 24, 'FontName', 'Times New Roman'); % Adjust the y-value as needed
end

% Create a legend for the first plot of each group
legend([h1, h2], {'Group A', 'Group B'}, 'FontName', 'Times New Roman', 'FontSize', 24);

% Set y-axis label
xlabel('Video Number', 'FontName', 'Times New Roman', 'FontSize', 24);
ylabel('View Duration (s)', 'FontName', 'Times New Roman', 'FontSize', 24);

dat_GA_view = [mod_1_data_GA(:,2); mod_2_data_GA(:,2); mod_3_data_GA(:,2); mod_4_data_GA(:,2)];
dat_GB_view = [mod_1_data_GB(:,2); mod_2_data_GB(:,2); mod_3_data_GB(:,2); mod_4_data_GB(:,2)];

dat_GA_watch = [mod_1_data_GA(:,4); mod_2_data_GA(:,4); mod_3_data_GA(:,4); mod_4_data_GA(:,4)];
dat_GB_watch = [mod_1_data_GB(:,4); mod_2_data_GB(:,4); mod_3_data_GB(:,4); mod_4_data_GB(:,4)];

dat_GA_view_buffer = [mod_1_data_GA_buffer(:,2); mod_2_data_GA_buffer(:,2); mod_3_data_GA_buffer(:,2); mod_4_data_GA_buffer(:,2)];
dat_GB_view_buffer = [mod_1_data_GB_buffer(:,2); mod_2_data_GB_buffer(:,2); mod_3_data_GB_buffer(:,2); mod_4_data_GB_buffer(:,2)];

dat_GA_watch_buffer = [mod_1_data_GA_buffer(:,4); mod_2_data_GA_buffer(:,4); mod_3_data_GA_buffer(:,4); mod_4_data_GA_buffer(:,4)];
dat_GB_watch_buffer = [mod_1_data_GB_buffer(:,4); mod_2_data_GB_buffer(:,4); mod_3_data_GB_buffer(:,4); mod_4_data_GB_buffer(:,4)];

dat_GA_view_final = [mod_1_data_GA_Final(:,2),; mod_2_data_GA_Final(:,2); mod_3_data_GA_Final(:,2); mod_4_data_GA_Final(:,2)];
dat_GB_view_final = [mod_1_data_GB_Final(:,2),; mod_2_data_GB_Final(:,2); mod_3_data_GB_Final(:,2); mod_4_data_GB_Final(:,2)];

dat_GA_watch_final = [mod_1_data_GA_Final(:,4),; mod_2_data_GA_Final(:,4); mod_3_data_GA_Final(:,4); mod_4_data_GA_Final(:,4)];
dat_GB_watch_final = [mod_1_data_GB_Final(:,4),; mod_2_data_GB_Final(:,4); mod_3_data_GB_Final(:,4); mod_4_data_GB_Final(:,4)];

[p_view, h_view, ~] = ranksum(dat_GA_view, dat_GB_view, 'tail', 'right', 'alpha', 0.01);
[p_watch, h_watch, ~] = ranksum(dat_GA_watch, dat_GB_watch, 'tail', 'both', 'alpha', 0.01);

[p_view_buffer, h_view_buffer, ~] = ranksum(dat_GA_view_buffer, dat_GB_view_buffer, 'tail', 'both', 'alpha', 0.01);
[p_watch_buffer, h_watch_buffer, ~] = ranksum(dat_GA_watch_buffer, dat_GB_watch_buffer, 'tail', 'both', 'alpha', 0.01);

[p_EC_view, h_EC_view, ~] = ranksum(dat_GA_view_final, dat_GB_view_final, 'tail', 'left', 'alpha', 0.01);
[p_EC_watch, h_EC_watch, ~] = ranksum(dat_GA_watch_final, dat_GB_watch_final, 'tail', 'left', 'alpha', 0.01);

% Store results in a cell array
results = {'Test', 'p-value', 'Reject Null';...
           'View', p_view, ternary(h_view, 'Yes', 'No');...
           'Watch', p_watch, ternary(h_watch, 'Yes', 'No');...
           'View Buffer', p_view_buffer, ternary(h_view_buffer, 'Yes', 'No');...
           'Watch Buffer', p_watch_buffer, ternary(h_watch_buffer, 'Yes', 'No');...
           'Final View', p_EC_view, ternary(h_EC_view, 'Yes', 'No');...
           'Final Watch', p_EC_watch, ternary(h_EC_watch, 'Yes', 'No')};

% Print the table
fprintf('Wilcoxon Rank Sum Test Results alpha = 0.01:\n');
fprintf('%-15s %-15s %-15s\n', results{1,:});
for i = 2:size(results,1)
    fprintf('%-15s %-15.4f %-15s\n', results{i,:});
end

first_vid_avg_GA = mean([mod_1_data_GA_buffer(1,1), mod_2_data_GA_buffer(1,1), mod_3_data_GA_buffer(1,1), mod_4_data_GA_buffer(1,1)]);
first_vid_avg_GB = mean([mod_1_data_GB_buffer(1,1), mod_2_data_GB_buffer(1,1), mod_3_data_GB_buffer(1,1), mod_4_data_GB_buffer(1,1)]);


avg_mods_GA = mean([mod_1_data_GA_buffer(2:end,1); mod_2_data_GA_buffer(2:end,1); mod_3_data_GA_buffer(2:(end-1),1); mod_4_data_GA_buffer(2:end,1)]);
avg_mods_GB = mean([mod_1_data_GB_buffer(2:end,1); mod_2_data_GB_buffer(2:end,1); mod_3_data_GB_buffer(2:(end-1),1); mod_4_data_GB_buffer(2:end,1)]);

end_vid_avg_GA = mean([mod_1_data_GA_buffer(end,1), mod_2_data_GA_buffer(end,1), mod_3_data_GA_buffer(end,1), mod_4_data_GA_buffer(end,1)]);
end_vid_avg_GB = mean([mod_1_data_GB_buffer(end,1), mod_2_data_GB_buffer(end,1), mod_3_data_GB_buffer(end,1), mod_4_data_GB_buffer(end,1)]);