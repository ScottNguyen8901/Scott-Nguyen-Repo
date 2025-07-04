function createBarCharts(mean_GA, mean_GB, p_values_GA, p_values_GB, alpha)
    % Create a new figure
    figure;

    % Define modules (1 to 4)
    modules = 1:4;
    hold on;

    % Set font properties
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);

    % Create bar charts for Group A and Group B
    bar(modules - 0.2, mean_GA, 0.4, 'b'); % Group A bars
    bar(modules + 0.2, mean_GB, 0.4, 'r'); % Group B bars

    % Add values and p-values for each bar
    for i = 1:length(modules)
        % Display bar values on top
        text(modules(i) - 0.2, mean_GA(i) + 0.1, sprintf('%.2f', mean_GA(i)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 16, 'Color', 'b', 'Interpreter', 'tex');
        text(modules(i) + 0.2, mean_GB(i) + 0.1, sprintf('%.2f', mean_GB(i)), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 16, 'Color', 'r', 'Interpreter', 'tex');

        % Display p-values and hypothesis decision
        text(modules(i) - 0.25, -2, sprintf('%.1e\n%s', p_values_GA(i), decision(p_values_GA(i), alpha)), ...
             'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', 'b');
        text(modules(i) + 0.25, -2, sprintf('%.1e\n%s', p_values_GB(i), decision(p_values_GB(i), alpha)), ...
             'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', 'r');
    end

    % Add labels and legend
    text(0.25, -1.5, 'P-Val:', 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', 'k');
    text(0.25, -3, 'Rej H_0:', 'HorizontalAlignment', 'center', 'FontSize', 16, 'Color', 'k');
    ylabel('Change in Score (%)', 'FontName', 'Times New Roman', 'FontSize', 24);
    legend('Group A', 'Group B', 'FontName', 'Times New Roman', 'FontSize', 24);

    % Set axis limits, ticks, and grid
    xlim([-0.25 5]);
    ylim([-5 33]);
    xticks(modules);
    xticklabels({'Module 1', 'Module 2', 'Module 3', 'Module 4'});
    set(gca, 'XTickLabel', get(gca, 'XTickLabel'), 'FontName', 'Times New Roman', 'FontSize', 24);
    grid on;
    axis square;

    hold off;
end
