% fileparser.m
function terrain_data = fileparser(filename)
    %FILEPARSER Reads terrain data from a .04 file using dlmread.
    %   terrain_data = FILEPARSER(filename) reads the specified .04 file and
    %   returns the terrain data as a matrix.

    % Check if the file exists
    if ~isfile(filename)
        error('File not found: %s', filename);
    end

    % Read the file contents using dlmread
    % Assuming the .04 file contains numeric data in two columns: distance and elevation
    terrain_data = dlmread(filename);
end
