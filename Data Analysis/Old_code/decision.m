function dec = decision(p_value, alpha)
    if p_value < alpha
        dec = 'Yes';
    else
        dec = 'No';
    end
end