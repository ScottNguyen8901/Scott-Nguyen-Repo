function date_vec = tle_epoch_2_datetime(tle_epoch)
    %
    % DESCRIPTION
    %   Convert TLE epoch time to a datetime vector.
    %
    % INPUTS        size    Type
    %   tle_epoch   (1,1)   Double  TLE epoch time in Julian Date format.
    %
    % OUTPUTS       size    Type
    %   datevec     (1,6)   Double  Datetime vector [year, month, day, hour, minute, second].
    %
    % NOTES
    %
    % FUNCTION

    % Extract the year and day of the year from the TLE epoch
    ymd = floor(tle_epoch);
    yr = fix(ymd / 1000);
    dofyr = mod(ymd, 1000);

    % Determine the correct year based on the value of yr
    if (yr < 57)
        year = yr + 2000;
    else
        year = yr + 1900;
    end

    % Calculate the decimal day of the year
    decidy = round((tle_epoch - ymd) * 10^8) / 10^8;
    temp = decidy * 24;  % Convert day to hours
    hh = fix(temp);
    
    temp = (temp - hh) * 60;  % Convert remaining to minutes
    mm = fix(temp);
    
    ss = (temp - mm) * 60;  % Remaining seconds

    % Determine the number of days in each month
    nd = eomday(year, 1:12);
    temp = cumsum(nd);

    % Find the month corresponding to the day of the year
    month = find(temp >= dofyr, 1);
    temp = temp(month) - dofyr;
    date = nd(month) - temp;

    % Create the output date vector
    date_vec = [year, month, date, hh, mm, ss];
end
