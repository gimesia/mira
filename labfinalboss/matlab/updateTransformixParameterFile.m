function updateTransformixParameterFile(filePath, parameter, value)
    % Read the file contents
    fileContents = fileread(filePath);
    
    % Define the parameter line format with parentheses
    parameterLine = sprintf('(%s %s)', parameter, value);
    
    % Check if the parameter already exists
    if contains(fileContents, parameter)
        % Use a regex to update the parameter line
        updatedContents = regexprep(fileContents, ['\(' parameter ' .*?\)'], parameterLine);
    else
        % Append the parameter at the end of the file with proper format
        updatedContents = [fileContents newline parameterLine];
    end
    
    % Write the updated contents back to the file
    fid = fopen(filePath, 'w');
    fwrite(fid, updatedContents);
    fclose(fid);
end
