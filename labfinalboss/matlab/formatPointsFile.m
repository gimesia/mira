function formatPointsFile(inputFile, outputFile)
    % Read the data from the input file
    data = readmatrix(inputFile);

    % Ensure the file contains exactly 300 rows
    numRows = size(data, 1);
    if numRows ~= 300
        error('The input file must contain exactly 300 rows of data.');
    end

    % Open the output file for writing
    fid = fopen(outputFile, 'w');
    if fid == -1
        error('Could not open the output file for writing.');
    end

    % Write the required headers
    fprintf(fid, 'index\n');
    fprintf(fid, '%d\n', numRows);

    % Write the data rows
    for i = 1:numRows
        fprintf(fid, '%d\t%d\t%d\n', data(i, :));
    end

    % Close the file
    fclose(fid);

    % fprintf('Formatted file saved as: %s\n', outputFile);
end
