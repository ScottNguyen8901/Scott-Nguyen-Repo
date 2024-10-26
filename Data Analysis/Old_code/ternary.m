% Helper function to convert logical value to 'Yes' or 'No'
function out = ternary(condition, true_val, false_val)
    if condition
        out = true_val;
    else
        out = false_val;
    end
end