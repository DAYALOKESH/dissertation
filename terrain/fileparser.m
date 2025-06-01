function terrain_data = fileparser(filename, varargin)
    % Check if file exists
    if ~exist(filename, 'file')
        error('File not found: %s', filename);
    end

    % Parse optional maxDistance parameter
    maxDistance = inf;
    if nargin > 1 && ~isempty(varargin{1})
        maxDistance = varargin{1};
    end
    
    try
        % Read terrain data using dlmread
        full_data = dlmread(filename);
        
        % Return all data or filter by distance
        if isinf(maxDistance)
            terrain_data = full_data;
        else
            % Get points up to maxDistance
            indices = full_data(:,1) <= maxDistance;
            terrain_data = full_data(indices, :);
            
            % Add one more point if available (for interpolation)
            if sum(indices) < size(full_data, 1)
                next_idx = find(indices, 1, 'last') + 1;
                if next_idx <= size(full_data, 1)
                    terrain_data = [terrain_data; full_data(next_idx, :)];
                end
            end
        end
    catch
        error('Error reading file: %s', filename);
    end
end