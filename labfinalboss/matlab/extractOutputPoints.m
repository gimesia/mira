function outputMatrix = extractOutputPoints(filename)
    % Open the file for reading
    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open the file: %s', filename);
    end

    % Initialize an empty array to hold the extracted points
    outputMatrix = [];

    % Read the file line by line
    while ~feof(fid)
        line = fgetl(fid);
        % Check if the line contains "OutputPoint"
        if contains(line, 'OutputIndexFixed =')
            % Extract the numbers after "OutputPoint = ["
            tokens = regexp(line, 'OutputIndexFixed = \[([-\d.e ]+)\]', 'tokens');
            if ~isempty(tokens)
                % Convert the string of numbers into a row vector
                point = str2num(tokens{1}{1}); %#ok<ST2NM>
                outputMatrix = [outputMatrix; point]; %#ok<AGROW>
            end
        end
    end

    % Close the file
    fclose(fid);

    % Ensure the output matrix has 300 rows
    if size(outputMatrix, 1) ~= 300
        error('The number of extracted points (%d) does not match the expected count (300).', size(outputMatrix, 1));
    end
end
