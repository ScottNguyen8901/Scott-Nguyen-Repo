function parse_ephemeris_data(folderPath, outputDirectory)
    % Initialize variables
    julianDate = [];
    x = [];
    y = [];
    z = [];
    vx = [];
    vy = [];
    vz = [];
    
    % Get list of .txt files in the folder matching the pattern 'planet_horizons_results.txt'
    fileList = dir(fullfile(folderPath, '*_horizons_results.txt'));
    totalFiles = length(fileList);  % Total number of files
    
    % Loop through each file in the folder
    for fileIdx = 1:totalFiles
        filePath = fullfile(folderPath, fileList(fileIdx).name);
        
        % Open and read the file
        fileID = fopen(filePath, 'r');
        data = fread(fileID, '*char')';
        fclose(fileID);
        
        % Extract the target body name from the file contents
        targetNameLine = regexp(data, 'Target body name: (.*?)(\r|\n)', 'tokens', 'once');
        if isempty(targetNameLine)
            disp(['Skipping file: ' fileList(fileIdx).name ' (target body name not found)']);
            continue;
        end
        
        targetName = strtrim(targetNameLine{1});
        
        % Normalize the target name by removing source info and any extra characters
        targetName = regexprep(targetName, '\{.*\}', ''); % Remove anything within { ... }
        targetName = regexprep(targetName, '[\(\[].*?[\)\]]', ''); % Remove parentheses or square brackets
        targetName = strrep(targetName, ' ', ''); % Remove spaces
        
        % Get the base name of the file (without the directory and extension)
        [~, fileBaseName, ~] = fileparts(fileList(fileIdx).name);
        
        % Normalize the file name by removing any underscores or spaces
        fileBaseName = strrep(fileBaseName, '_', '');  % Remove underscores from the file name
        fileBaseName = strrep(fileBaseName, ' ', '');  % Remove spaces
        
        if contains(fileBaseName, targetName, 'IgnoreCase', true)
            disp(['File Name matches the Target Name: ', targetName]);
        else
            disp(['Skipping file: ' fileList(fileIdx).name ' (file name does not match target body name)']);
            continue;
        end

        % Split data into lines
        dataLines = splitlines(data);
        
        % Parse data lines
        for i = 1:length(dataLines)
            line = strtrim(dataLines{i});
            if contains(line, '=')
                % Extract Julian Date (numeric part before the '=' symbol)
                jdStr = extractBefore(line, ' =');
                jd = str2double(strtrim(jdStr));
                if ~isnan(jd)
                    julianDate(end + 1, 1) = jd;
                end
                
                % Parse position and velocity components
                if contains(line, 'X =')
                    xVal = str2double(extractBetween(line, 'X =', 'Y ='));
                    yVal = str2double(extractBetween(line, 'Y =', 'Z ='));
                    zVal = str2double(extractAfter(line, 'Z ='));
                    x(end + 1, 1) = xVal;
                    y(end + 1, 1) = yVal;
                    z(end + 1, 1) = zVal;
                elseif contains(line, 'VX=') 
                    vxVal = str2double(extractBetween(line, 'VX=', 'VY='));
                    vyVal = str2double(extractBetween(line, 'VY=', 'VZ='));
                    vzVal = str2double(extractAfter(line, 'VZ='));
                    vx(end + 1, 1) = vxVal;
                    vy(end + 1, 1) = vyVal;
                    vz(end + 1, 1) = vzVal;
                end
            end
        end
        
        % Create table with all extracted data
        dataTable = table(julianDate, x, y, z, vx, vy, vz, ...
                          'VariableNames', {'JulianDate', 'X', 'Y', 'Z', 'VX', 'VY', 'VZ'});
        
        % Define output file path based on target name
        outputFileName = strcat(targetName, '_ephemeris_data.txt');
        outputFilePath = fullfile(outputDirectory, outputFileName);
        
        % Write the table to a .txt file
        writetable(dataTable, outputFilePath, 'Delimiter', '\t');
        
        % Calculate the progress percentage
        progress = (fileIdx / totalFiles) * 100;
        fprintf('Processing: %.2f%%\n', progress);
    end
end