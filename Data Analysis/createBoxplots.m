function createBoxplots(dataA, dataB, pos_offset, colorA, colorB)
    figure; % Create a new figure
    
    % Define positions for the boxplots
    pos_A = (1:4) - pos_offset;
    pos_B = pos_A + 2 * pos_offset;

    % First subplot (PC Scores)
    subplot1 = axes('Position', [0.1, 0.1, 0.35, 0.8]); % Position the first subplot
    h_A = boxplot(dataA, 'Labels', {'Module 1', 'Module 2', 'Module 3', 'Module 4'}, ...
                  'Colors', colorA, 'Widths', 0.3, 'Positions', pos_A);
    hold on;

    % Plot Group B Boxplot for PC
    h_B = boxplot(dataB, 'Labels', {}, 'Colors', colorB, 'Widths', 0.3, 'Positions', pos_B);
    
    % Set axis properties for PC plot
    ylabel('Scores (%)', 'FontName', 'Times New Roman', 'FontSize', 24);
    set(gca, 'XTick', 1:4, 'XTickLabel', {'Module 1', 'Module 2', 'Module 3', 'Module 4'}, ...
             'FontSize', 24, 'FontName', 'Times New Roman');
    grid on; axis square;
    set(findobj(h_A, 'type', 'line'), 'LineWidth', 2);
    set(findobj(h_B, 'type', 'line'), 'LineWidth', 2);
    set(findobj(h_A, 'Tag', 'Outliers'), 'MarkerEdgeColor', colorA);
    set(findobj(h_B, 'Tag', 'Outliers'), 'MarkerEdgeColor', colorB);

    % Add mean values as text for PC plot
    for i = 1:4
        text(pos_A(i), mean(dataA(:, i)) + 5, sprintf('%.2f', mean(dataA(:, i))), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'FontSize', 10, 'Color', colorA, 'FontName', 'Times New Roman');
        text(pos_B(i), mean(dataB(:, i)) + 5, sprintf('%.2f', mean(dataB(:, i))), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'FontSize', 10, 'Color', colorB, 'FontName', 'Times New Roman');
    end
    
    % Second subplot (MC Scores)
    subplot2 = axes('Position', [0.55, 0.1, 0.35, 0.8]); % Position the second subplot
    h_A = boxplot(dataA, 'Labels', {'Module 1', 'Module 2', 'Module 3', 'Module 4'}, ...
                  'Colors', colorA, 'Widths', 0.3, 'Positions', pos_A);
    hold on;

    % Plot Group B Boxplot for MC
    h_B = boxplot(dataB, 'Labels', {}, 'Colors', colorB, 'Widths', 0.3, 'Positions', pos_B);
    
    % Set axis properties for MC plot
    ylabel('Scores (%)', 'FontName', 'Times New Roman', 'FontSize', 24);
    set(gca, 'XTick', 1:4, 'XTickLabel', {'Module 1', 'Module 2', 'Module 3', 'Module 4'}, ...
             'FontSize', 24, 'FontName', 'Times New Roman');
    grid on; axis square;
    set(findobj(h_A, 'type', 'line'), 'LineWidth', 2);
    set(findobj(h_B, 'type', 'line'), 'LineWidth', 2);
    set(findobj(h_A, 'Tag', 'Outliers'), 'MarkerEdgeColor', colorA);
    set(findobj(h_B, 'Tag', 'Outliers'), 'MarkerEdgeColor', colorB);

    % Add mean values as text for MC plot
    for i = 1:4
        text(pos_A(i), mean(dataA(:, i)) + 5, sprintf('%.2f', mean(dataA(:, i))), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'FontSize', 10, 'Color', colorA, 'FontName', 'Times New Roman');
        text(pos_B(i), mean(dataB(:, i)) + 5, sprintf('%.2f', mean(dataB(:, i))), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'FontSize', 10, 'Color', colorB, 'FontName', 'Times New Roman');
    end

    % Add legend for both subplots
    legend([h_A(5), h_B(5)], 'Group A', 'Group B', 'FontName', 'Times New Roman', 'FontSize', 24, 'Location', 'northeast');
end
