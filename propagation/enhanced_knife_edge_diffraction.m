function [diff_loss, field_with_diff, diffraction_info] = enhanced_knife_edge_diffraction(varargin)
%ENHANCED_KNIFE_EDGE_DIFFRACTION Calculate diffraction loss using advanced knife edge models
%
%   [diff_loss, field_with_diff, diffraction_info] = ENHANCED_KNIFE_EDGE_DIFFRACTION()
%   [diff_loss, field_with_diff, diffraction_info] = ENHANCED_KNIFE_EDGE_DIFFRACTION('param', value, ...)
%
%   This function implements advanced knife edge diffraction models to calculate
%   the additional loss due to terrain obstacles in the propagation path.
%   It supports single and multiple knife edges and implements several diffraction methods.
%
%   Parameters:
%       'frequency'     - Frequency in MHz (default: 970)
%       'txHeight'      - Transmitter height in meters (default: 50)
%       'rxHeight'      - Receiver height in meters (default: 1.5)
%       'txPosition'    - Transmitter position [x,y] (default: [0,50])
%       'rxPosition'    - Receiver position [x,y] (default: [700,1.5])
%       'terrain'       - Terrain data as [x,y] array (required)
%       'fieldStrength' - Original field strength in dB (required)
%       'method'        - Diffraction method: 'itu', 'deygout', 'epstein-peterson', 'bullington' (default: 'itu')
%       'maxEdges'      - Maximum number of knife edges to consider (default: 3)
%       'outputDir'     - Directory to save output files (default: './results')
%       'plotResults'   - Whether to plot results (default: true)
%
%   Returns:
%       diff_loss       - Diffraction loss in dB
%       field_with_diff - Field strength with diffraction in dB
%       diffraction_info- Structure with detailed diffraction information
%
%   Current Date and Time (UTC - YYYY-MM-DD HH:MM:SS formatted): 2025-06-01 15:54:12
%   Current User's Login: DAYALOKESH

    % Parse input parameters
    p = inputParser;
    p.addParameter('frequency', 970, @isnumeric);    % MHz
    p.addParameter('txHeight', 50, @isnumeric);      % m
    p.addParameter('rxHeight', 1.5, @isnumeric);     % m
    p.addParameter('txPosition', [0, 50], @isnumeric);  % [x,y]
    p.addParameter('rxPosition', [700, 1.5], @isnumeric); % [x,y]
    p.addParameter('terrain', [], @isnumeric);       % Terrain data
    p.addParameter('fieldStrength', [], @isnumeric); % Original field strength
    p.addParameter('method', 'itu', @ischar);        % Diffraction method
    p.addParameter('maxEdges', 3, @isnumeric);       % Max number of knife edges
    p.addParameter('outputDir', './results', @ischar); % Output directory
    p.addParameter('plotResults', true, @islogical); % Plot results flag
    p.parse(varargin{:});
    
    params = p.Results;
    
    % Check required parameters
    if isempty(params.terrain) || isempty(params.fieldStrength)
        error('Terrain data and field strength are required parameters');
    end
    
    % Ensure output directory exists
    if ~exist(params.outputDir, 'dir')
        mkdir(params.outputDir);
        fprintf('Created output directory: %s\n', params.outputDir);
    end
    
    % Get terrain data
    terrain_x = params.terrain(:,1);
    terrain_y = params.terrain(:,2);
    
    % Get transmitter and receiver positions
    tx_x = params.txPosition(1);
    tx_y = params.txPosition(2);
    rx_x = params.rxPosition(1);
    rx_y = params.rxPosition(2);
    
    % Calculate wavelength
    lambda = 299.792458 / params.frequency;  % Wavelength in meters
    
    % Number of points to analyze
    num_points = length(terrain_x);
    
    fprintf('\n=== ENHANCED KNIFE EDGE DIFFRACTION ANALYSIS ===\n');
    fprintf('Frequency: %.2f MHz (Wavelength: %.4f m)\n', params.frequency, lambda);
    fprintf('TX Position: (%.2f, %.2f) m\n', tx_x, tx_y);
    fprintf('RX Position: (%.2f, %.2f) m\n', rx_x, rx_y);
    fprintf('Using diffraction method: %s\n', upper(params.method));
    fprintf('Analyzing %d terrain points for up to %d knife edges\n', num_points, params.maxEdges);
    
    % Calculate the direct line-of-sight path between TX and RX
    % Use parametric equation of a line to handle vertical paths
    tx_to_rx_dist = sqrt((rx_x - tx_x)^2 + (rx_y - tx_y)^2);
    num_los_points = 1000; % High resolution for LOS analysis
    
    % Generate points along the line from TX to RX
    t = linspace(0, 1, num_los_points);
    los_x = tx_x + (rx_x - tx_x) * t;
    los_y = tx_y + (rx_y - tx_y) * t;
    
    % Interpolate terrain heights at the LOS x-coordinates
    terrain_interp_y = interp1(terrain_x, terrain_y, los_x, 'linear', 'extrap');
    
    % Calculate clearance (height above terrain) for each point along LOS
    clearance = los_y - terrain_interp_y;
    
    % Identify potential knife edges (obstacles) where clearance is negative or minimal
    potential_edges = [];
    min_clearance_threshold = -lambda; % Consider points with negative clearance as potential edges
    
    for i = 2:(num_los_points-1)
        % Check if this point is a local minimum in clearance
        if (clearance(i) <= clearance(i-1) && clearance(i) <= clearance(i+1)) || clearance(i) < min_clearance_threshold
            % Calculate the angle from horizontal for the previous, current, and next point
            angle_prev = atan2(los_y(i) - los_y(i-1), los_x(i) - los_x(i-1)) * 180/pi;
            angle_next = atan2(los_y(i+1) - los_y(i), los_x(i+1) - los_x(i)) * 180/pi;
            
            % If there's a significant angle change, it's a potential edge
            if abs(angle_next - angle_prev) > 1 || clearance(i) < min_clearance_threshold
                % Store x, y, clearance, and index of potential edge
                potential_edges = [potential_edges; los_x(i), terrain_interp_y(i), clearance(i), i];
            end
        end
    end
    
    % If no potential edges found, add the point with minimum clearance
    if isempty(potential_edges)
        [min_clear, min_idx] = min(clearance);
        potential_edges = [los_x(min_idx), terrain_interp_y(min_idx), min_clear, min_idx];
    end
    
    % Sort edges by clearance (most significant obstruction first)
    [~, sort_idx] = sort(potential_edges(:,3));
    potential_edges = potential_edges(sort_idx,:);
    
    % Limit to maximum number of edges to consider
    num_edges = min(size(potential_edges, 1), params.maxEdges);
    potential_edges = potential_edges(1:num_edges, :);
    
    fprintf('Identified %d significant knife edges\n', num_edges);
    
    % Calculate diffraction loss based on the selected method
    switch lower(params.method)
        case 'itu'
            % ITU-R P.526 method
            [diff_loss, edge_losses] = calculate_itu_diffraction(potential_edges, tx_x, tx_y, rx_x, rx_y, lambda);
            
        case 'deygout'
            % Deygout method for multiple knife edges
            [diff_loss, edge_losses] = calculate_deygout_diffraction(potential_edges, tx_x, tx_y, rx_x, rx_y, lambda);
            
        case 'epstein-peterson'
            % Epstein-Peterson method for multiple knife edges
            [diff_loss, edge_losses] = calculate_epstein_peterson_diffraction(potential_edges, tx_x, tx_y, rx_x, rx_y, lambda);
            
        case 'bullington'
            % Bullington method (equivalent single knife edge)
            [diff_loss, edge_losses] = calculate_bullington_diffraction(potential_edges, tx_x, tx_y, rx_x, rx_y, lambda);
            
        otherwise
            % Default to ITU-R P.526
            fprintf('Unknown method "%s", defaulting to ITU-R P.526\n', params.method);
            [diff_loss, edge_losses] = calculate_itu_diffraction(potential_edges, tx_x, tx_y, rx_x, rx_y, lambda);
    end
    
    % Apply diffraction loss to the original field strength
    field_with_diff = params.fieldStrength - diff_loss;
    
    % Display results
    fprintf('\nDiffraction Analysis Results:\n');
    fprintf('- Total diffraction loss: %.2f dB\n', diff_loss);
    fprintf('- Number of edges analyzed: %d\n', num_edges);
    
    for i = 1:num_edges
        fprintf('- Edge %d at (%.1f, %.1f) m: %.2f dB loss, clearance: %.2f m\n', ...
                i, potential_edges(i,1), potential_edges(i,2), edge_losses(i), potential_edges(i,3));
    end
    
    % Prepare return structure with detailed diffraction information
    diffraction_info = struct();
    diffraction_info.method = params.method;
    diffraction_info.total_loss = diff_loss;
    diffraction_info.edges = potential_edges;
    diffraction_info.edge_losses = edge_losses;
    diffraction_info.los_x = los_x;
    diffraction_info.los_y = los_y;
    diffraction_info.clearance = clearance;
    diffraction_info.terrain_interp_y = terrain_interp_y;
    diffraction_info.lambda = lambda;
    
    % Save results to file
    outputFile = fullfile(params.outputDir, ['KnifeEdge_', upper(params.method), '_Results.txt']);
    saveResults(params, diffraction_info, outputFile, field_with_diff);
    
    % Plot results if requested
    if params.plotResults
        plotDiffractionResults(params, diffraction_info, field_with_diff, outputFile);
    end
end

function [diff_loss, edge_losses] = calculate_itu_diffraction(edges, tx_x, tx_y, rx_x, rx_y, lambda)
    % Calculate diffraction loss using ITU-R P.526 method
    num_edges = size(edges, 1);
    edge_losses = zeros(num_edges, 1);
    
    % For each knife edge, calculate Fresnel-Kirchhoff diffraction parameter v
    for i = 1:num_edges
        edge_x = edges(i, 1);
        edge_y = edges(i, 2);
        
        % Calculate direct distance from TX to RX
        d_direct = sqrt((rx_x - tx_x)^2 + (rx_y - tx_y)^2);
        
        % Calculate distances from TX to edge and edge to RX
        d1 = sqrt((edge_x - tx_x)^2 + (edge_y - tx_y)^2);
        d2 = sqrt((rx_x - edge_x)^2 + (rx_y - edge_y)^2);
        
        % Calculate the height of the direct path at the edge position
        % Using similar triangles
        h_direct = tx_y + (rx_y - tx_y) * (edge_x - tx_x) / (rx_x - tx_x);
        
        % Calculate the height difference
        h_diff = edge_y - h_direct;
        
        % Calculate Fresnel-Kirchhoff diffraction parameter v
        v = h_diff * sqrt(2 * d1 * d2 / (lambda * d_direct * (d1 + d2)));
        
        % Calculate diffraction loss using ITU-R P.526 formula
        if v > -0.78
            % For v > -0.78, use the approximation formula
            J_v = 6.9 + 20 * log10(sqrt((v - 0.1)^2 + 1) + v - 0.1);
        else
            % For v ≤ -0.78, no diffraction loss
            J_v = 0;
        end
        
        edge_losses(i) = J_v;
    end
    
    % Total diffraction loss is the sum of individual edge losses
    diff_loss = sum(edge_losses);
end

function [diff_loss, edge_losses] = calculate_deygout_diffraction(edges, tx_x, tx_y, rx_x, rx_y, lambda)
    % Calculate diffraction loss using Deygout method (main edge + subsidiary edges)
    num_edges = size(edges, 1);
    edge_losses = zeros(num_edges, 1);
    
    if num_edges == 1
        % Single edge case - use ITU method
        [diff_loss, edge_losses] = calculate_itu_diffraction(edges, tx_x, tx_y, rx_x, rx_y, lambda);
        return;
    end
    
    % Sort edges by x-coordinate
    [~, sort_idx] = sort(edges(:, 1));
    edges = edges(sort_idx, :);
    
    % Find the main edge (maximum v parameter)
    v_values = zeros(num_edges, 1);
    for i = 1:num_edges
        edge_x = edges(i, 1);
        edge_y = edges(i, 2);
        
        % Calculate direct distance from TX to RX
        d_direct = sqrt((rx_x - tx_x)^2 + (rx_y - tx_y)^2);
        
        % Calculate distances from TX to edge and edge to RX
        d1 = sqrt((edge_x - tx_x)^2 + (edge_y - tx_y)^2);
        d2 = sqrt((rx_x - edge_x)^2 + (rx_y - edge_y)^2);
        
        % Calculate the height of the direct path at the edge position
        h_direct = tx_y + (rx_y - tx_y) * (edge_x - tx_x) / (rx_x - tx_x);
        
        % Calculate the height difference
        h_diff = edge_y - h_direct;
        
        % Calculate Fresnel-Kirchhoff diffraction parameter v
        v_values(i) = h_diff * sqrt(2 * d1 * d2 / (lambda * d_direct * (d1 + d2)));
    end
    
    % Find the main edge (maximum v parameter)
    [~, main_idx] = max(v_values);
    main_edge = edges(main_idx, :);
    main_x = main_edge(1);
    main_y = main_edge(2);
    
    % Calculate main edge diffraction using ITU method
    main_edge_loss = calculate_single_edge_loss(main_edge, tx_x, tx_y, rx_x, rx_y, lambda);
    edge_losses(main_idx) = main_edge_loss;
    
    % Calculate subsidiary edges diffraction
    % Edges between TX and main edge
    if main_idx > 1
        sub_edges_tx = edges(1:main_idx-1, :);
        [sub_loss_tx, sub_losses_tx] = calculate_itu_diffraction(sub_edges_tx, tx_x, tx_y, main_x, main_y, lambda);
        edge_losses(1:main_idx-1) = sub_losses_tx;
    else
        sub_loss_tx = 0;
    end
    
    % Edges between main edge and RX
    if main_idx < num_edges
        sub_edges_rx = edges(main_idx+1:end, :);
        [sub_loss_rx, sub_losses_rx] = calculate_itu_diffraction(sub_edges_rx, main_x, main_y, rx_x, rx_y, lambda);
        edge_losses(main_idx+1:end) = sub_losses_rx;
    else
        sub_loss_rx = 0;
    end
    
    % Calculate total diffraction loss using Deygout method
    diff_loss = main_edge_loss + sub_loss_tx + sub_loss_rx;
    
    % Apply correction factor to avoid overestimation
    if num_edges > 2
        % Apply correction based on ITU-R P.526-15
        correction = 10 * log10(0.5);  % -3 dB correction
        diff_loss = diff_loss + correction;
    end
end

function [diff_loss, edge_losses] = calculate_epstein_peterson_diffraction(edges, tx_x, tx_y, rx_x, rx_y, lambda)
    % Calculate diffraction loss using Epstein-Peterson method
    num_edges = size(edges, 1);
    edge_losses = zeros(num_edges, 1);
    
    if num_edges == 1
        % Single edge case - use ITU method
        [diff_loss, edge_losses] = calculate_itu_diffraction(edges, tx_x, tx_y, rx_x, rx_y, lambda);
        return;
    end
    
    % Sort edges by x-coordinate
    [~, sort_idx] = sort(edges(:, 1));
    edges = edges(sort_idx, :);
    
    % Calculate diffraction for each edge considering only adjacent edges
    total_loss = 0;
    
    for i = 1:num_edges
        edge_x = edges(i, 1);
        edge_y = edges(i, 2);
        
        if i == 1
            % First edge: TX to edge to second edge
            next_x = edges(i+1, 1);
            next_y = edges(i+1, 2);
            edge_loss = calculate_single_edge_loss([edge_x, edge_y], tx_x, tx_y, next_x, next_y, lambda);
        elseif i == num_edges
            % Last edge: previous edge to last edge to RX
            prev_x = edges(i-1, 1);
            prev_y = edges(i-1, 2);
            edge_loss = calculate_single_edge_loss([edge_x, edge_y], prev_x, prev_y, rx_x, rx_y, lambda);
        else
            % Middle edges: previous edge to current edge to next edge
            prev_x = edges(i-1, 1);
            prev_y = edges(i-1, 2);
            next_x = edges(i+1, 1);
            next_y = edges(i+1, 2);
            edge_loss = calculate_single_edge_loss([edge_x, edge_y], prev_x, prev_y, next_x, next_y, lambda);
        end
        
        edge_losses(i) = edge_loss;
        total_loss = total_loss + edge_loss;
    end
    
    diff_loss = total_loss;
end

function [diff_loss, edge_losses] = calculate_bullington_diffraction(edges, tx_x, tx_y, rx_x, rx_y, lambda)
    % Calculate diffraction loss using Bullington method (equivalent single knife edge)
    num_edges = size(edges, 1);
    
    % For single edge, just use ITU method
    if num_edges == 1
        [diff_loss, edge_losses] = calculate_itu_diffraction(edges, tx_x, tx_y, rx_x, rx_y, lambda);
        return;
    end
    
    % Sort edges by x-coordinate
    [~, sort_idx] = sort(edges(:, 1));
    edges = edges(sort_idx, :);
    
    % Find the virtual knife edge (intersection of TX and RX horizon rays)
    % Calculate slopes from TX to each edge
    slopes_tx = zeros(num_edges, 1);
    for i = 1:num_edges
        slopes_tx(i) = (edges(i,2) - tx_y) / (edges(i,1) - tx_x);
    end
    
    % Find maximum slope from TX
    [max_slope_tx, max_idx_tx] = max(slopes_tx);
    tx_horizon_edge = edges(max_idx_tx, :);
    
    % Calculate slopes from RX to each edge
    slopes_rx = zeros(num_edges, 1);
    for i = 1:num_edges
        slopes_rx(i) = (edges(i,2) - rx_y) / (edges(i,1) - rx_x);
    end
    
    % Find maximum slope from RX (considering direction)
    [max_slope_rx, max_idx_rx] = max(-slopes_rx);
    max_slope_rx = -max_slope_rx;
    rx_horizon_edge = edges(max_idx_rx, :);
    
    % Check if horizon rays intersect
    if max_slope_tx >= max_slope_rx
        % Rays intersect - calculate virtual knife edge position
        % y = mx + b equations for both lines
        b_tx = tx_y - max_slope_tx * tx_x;
        b_rx = rx_y - max_slope_rx * rx_x;
        
        % Intersection point
        x_virtual = (b_rx - b_tx) / (max_slope_tx - max_slope_rx);
        y_virtual = max_slope_tx * x_virtual + b_tx;
        
        % Create virtual edge
        virtual_edge = [x_virtual, y_virtual, 0, 0];
        
        % Calculate diffraction loss for virtual edge
        [diff_loss, edge_losses] = calculate_itu_diffraction(virtual_edge, tx_x, tx_y, rx_x, rx_y, lambda);
    else
        % Rays don't intersect - no diffraction loss
        diff_loss = 0;
        edge_losses = zeros(num_edges, 1);
    end
end

function edge_loss = calculate_single_edge_loss(edge, tx_x, tx_y, rx_x, rx_y, lambda)
    % Calculate diffraction loss for a single edge using ITU-R P.526 method
    edge_x = edge(1);
    edge_y = edge(2);
    
    % Calculate direct distance from TX to RX
    d_direct = sqrt((rx_x - tx_x)^2 + (rx_y - tx_y)^2);
    
    % Calculate distances from TX to edge and edge to RX
    d1 = sqrt((edge_x - tx_x)^2 + (edge_y - tx_y)^2);
    d2 = sqrt((rx_x - edge_x)^2 + (rx_y - edge_y)^2);
    
    % Calculate the height of the direct path at the edge position
    h_direct = tx_y + (rx_y - tx_y) * (edge_x - tx_x) / (rx_x - tx_x);
    
    % Calculate the height difference
    h_diff = edge_y - h_direct;
    
    % Calculate Fresnel-Kirchhoff diffraction parameter v
    v = h_diff * sqrt(2 * d1 * d2 / (lambda * d_direct * (d1 + d2)));
    
    % Calculate diffraction loss using ITU-R P.526 formula
    if v > -0.78
        % For v > -0.78, use the approximation formula
        edge_loss = 6.9 + 20 * log10(sqrt((v - 0.1)^2 + 1) + v - 0.1);
    else
        % For v ≤ -0.78, no diffraction loss
        edge_loss = 0;
    end
end

function saveResults(params, diffraction_info, outputFile, field_with_diff)
    % Save results to a file
    fid = fopen(outputFile, 'w');
    
    if fid == -1
        warning('Could not open output file for writing.');
        return;
    end
    
    % Write header with timestamp
    fprintf(fid, '# Enhanced Knife Edge Diffraction Results\n');
    fprintf(fid, '# Current Date and Time (UTC - YYYY-MM-DD HH:MM:SS formatted): 2025-06-01 15:54:12\n');
    fprintf(fid, '# Current User''s Login: DAYALOKESH\n\n');
    
    % Write parameters
    fprintf(fid, 'Parameters:\n');
    fprintf(fid, '- Frequency: %.2f MHz (Wavelength: %.4f m)\n', params.frequency, diffraction_info.lambda);
    fprintf(fid, '- Transmitter position: [%.2f, %.2f] m\n', params.txPosition(1), params.txPosition(2));
    fprintf(fid, '- Receiver position: [%.2f, %.2f] m\n', params.rxPosition(1), params.rxPosition(2));
    fprintf(fid, '- Diffraction method: %s\n', upper(params.method));
    fprintf(fid, '- Maximum edges considered: %d\n\n', params.maxEdges);
    
    % Write diffraction analysis results
    fprintf(fid, 'Diffraction Analysis Results:\n');
    fprintf(fid, '- Total diffraction loss: %.2f dB\n', diffraction_info.total_loss);
    fprintf(fid, '- Number of edges analyzed: %d\n\n', size(diffraction_info.edges, 1));
    
    % Write information about each knife edge
    fprintf(fid, 'Knife Edge Details:\n');
    for i = 1:size(diffraction_info.edges, 1)
        fprintf(fid, '- Edge %d at (%.2f, %.2f) m: %.2f dB loss, clearance: %.2f m\n', ...
                i, diffraction_info.edges(i,1), diffraction_info.edges(i,2), ...
                diffraction_info.edge_losses(i), diffraction_info.edges(i,3));
    end
    fprintf(fid, '\n');
    
    % Write detailed results for each point
    fprintf(fid, 'Detailed Results:\n');
    fprintf(fid, 'Distance(m)\tField_With_Diffraction(dB)\n');
    
    for i = 1:length(params.fieldStrength)
        % If terrain contains all points
        if i <= size(params.terrain, 1)
            fprintf(fid, '%.6f\t%.6f\n', params.terrain(i, 1), field_with_diff(i));
        end
    end
    
    % Close file
    fclose(fid);
    fprintf('Results saved to: %s\n', outputFile);
    
    % Also save data file for potential further analysis
    [filepath, name, ~] = fileparts(outputFile);
    matfile = fullfile(filepath, [name, '.mat']);
    save(matfile, 'params', 'diffraction_info', 'field_with_diff');
    fprintf('MATLAB data saved to: %s\n', matfile);
end

function plotDiffractionResults(params, diffraction_info, field_with_diff, outputFile)
    % Create figure for enhanced diffraction results
    figure('Name', ['Enhanced Knife Edge Diffraction Analysis - ', upper(params.method)], ...
           'Position', [100, 100, 1000, 800]);
    
    % Plot 1: Terrain profile with line-of-sight path and knife edges
    subplot(3, 1, 1);
    
    % Plot terrain
    if size(params.terrain, 1) <= 1000
        terrain_x = params.terrain(:,1);
        terrain_y = params.terrain(:,2);
    else
        % Downsample terrain for better plot performance
        downsample_factor = floor(size(params.terrain, 1) / 1000);
        terrain_x = params.terrain(1:downsample_factor:end,1);
        terrain_y = params.terrain(1:downsample_factor:end,2);
    end
    
    plot(terrain_x, terrain_y, 'k-', 'LineWidth', 1.5); hold on;
    
    % Plot line-of-sight path
    plot(diffraction_info.los_x, diffraction_info.los_y, 'b--', 'LineWidth', 1.5);
    
    % Plot TX and RX
    plot(params.txPosition(1), params.txPosition(2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    plot(params.rxPosition(1), params.rxPosition(2), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
    
    % Plot knife edges
    edge_colors = lines(size(diffraction_info.edges, 1));  % Get colormap for edges
    for i = 1:size(diffraction_info.edges, 1)
        edge_x = diffraction_info.edges(i,1);
        edge_y = diffraction_info.edges(i,2);
        plot(edge_x, edge_y, 'o', 'Color', edge_colors(i,:), 'MarkerSize', 10, ...
             'MarkerFaceColor', edge_colors(i,:));
        
        % Add edge number label
        text(edge_x + 5, edge_y + 5, sprintf('Edge %d\n%.1f dB', i, diffraction_info.edge_losses(i)), ...
             'FontWeight', 'bold', 'Color', edge_colors(i,:));
    end
    
    % Add 1st Fresnel zone ellipse visualization
    % Calculate ellipse for 1st Fresnel zone
    % Simplified visualization - just show radius at each knife edge
    for i = 1:size(diffraction_info.edges, 1)
        edge_x = diffraction_info.edges(i,1);
        edge_y = diffraction_info.edges(i,2);
        
        % Calculate direct distance from TX to RX
        d_direct = sqrt((params.rxPosition(1) - params.txPosition(1))^2 + ...
                         (params.rxPosition(2) - params.txPosition(2))^2);
        
        % Calculate distances from TX to edge and edge to RX
        d1 = sqrt((edge_x - params.txPosition(1))^2 + (edge_y - params.txPosition(2))^2);
        d2 = sqrt((params.rxPosition(1) - edge_x)^2 + (params.rxPosition(2) - edge_y)^2);
        
        % Calculate 1st Fresnel zone radius at this point
        F1 = sqrt(diffraction_info.lambda * d1 * d2 / d_direct);
        
        % Draw Fresnel zone radius as vertical line
        h_direct = params.txPosition(2) + (params.rxPosition(2) - params.txPosition(2)) * ...
                  (edge_x - params.txPosition(1)) / (params.rxPosition(1) - params.txPosition(1));
        
        % Draw vertical line for Fresnel zone radius (both above and below LOS)
        plot([edge_x, edge_x], [h_direct - F1, h_direct + F1], ':', 'Color', edge_colors(i,:), 'LineWidth', 1.5);
    end
    
    xlabel('Distance (m)', 'FontWeight', 'bold');
    ylabel('Height (m)', 'FontWeight', 'bold');
    title(sprintf('Terrain Profile with Line-of-Sight Path and Knife Edges (%s Method)', upper(params.method)), ...
          'FontWeight', 'bold');
    
    % Create legend
    legend_items = {'Terrain', 'Line-of-Sight', 'Transmitter', 'Receiver'};
    for i = 1:size(diffraction_info.edges, 1)
        legend_items{end+1} = sprintf('Edge %d', i);
    end
    legend(legend_items, 'Location', 'best');
    grid on;
    
    % Plot 2: Clearance profile
    subplot(3, 1, 2);
    
    % Calculate clearance in wavelengths
    clearance_wavelengths = diffraction_info.clearance / diffraction_info.lambda;
    
    % Plot clearance
    plot(diffraction_info.los_x, clearance_wavelengths, 'b-', 'LineWidth', 2); hold on;
    
    % Plot critical clearance thresholds
    plot(diffraction_info.los_x, zeros(size(diffraction_info.los_x)), 'k--');  % Zero clearance line
    plot(diffraction_info.los_x, 0.6 * ones(size(diffraction_info.los_x)), 'g--');  % 0.6λ line (60% Fresnel zone)
    plot(diffraction_info.los_x, -0.6 * ones(size(diffraction_info.los_x)), 'r--');  % -0.6λ line
    
    % Plot knife edge clearances
    for i = 1:size(diffraction_info.edges, 1)
        edge_x = diffraction_info.edges(i,1);
        edge_idx = find(abs(diffraction_info.los_x - edge_x) < 1, 1);
        if ~isempty(edge_idx)
            plot(edge_x, clearance_wavelengths(edge_idx), 'o', 'Color', edge_colors(i,:), ...
                 'MarkerSize', 10, 'MarkerFaceColor', edge_colors(i,:));
        end
    end
    
    xlabel('Distance (m)', 'FontWeight', 'bold');
    ylabel('Clearance (wavelengths)', 'FontWeight', 'bold');
    title('Clearance Above Line-of-Sight Path (in Wavelengths)', 'FontWeight', 'bold');
    legend('Clearance', 'Zero Clearance', '60% Fresnel Zone (Good)', 'Obstruction (-60%)', 'Location', 'best');
    grid on;
    
    % Plot 3: Field strength with diffraction
    subplot(3, 1, 3);
    
    % Downsample terrain if needed
    if length(params.fieldStrength) > 1000
        downsample_factor = floor(length(params.fieldStrength) / 1000);
        plot_indices = 1:downsample_factor:length(params.fieldStrength);
    else
        plot_indices = 1:length(params.fieldStrength);
    end
    
    % Get x-coordinates
    if size(params.terrain, 1) >= length(params.fieldStrength)
        plot_x = params.terrain(plot_indices, 1);
    else
        plot_x = linspace(params.terrain(1,1), params.terrain(end,1), length(params.fieldStrength));
        plot_x = plot_x(plot_indices);
    end
    
    % Plot field strength
    plot(plot_x, params.fieldStrength(plot_indices), 'b-', 'LineWidth', 2); hold on;
    plot(plot_x, field_with_diff(plot_indices), 'r-', 'LineWidth', 2);
    
    % Plot knife edge positions
    for i = 1:size(diffraction_info.edges, 1)
        edge_x = diffraction_info.edges(i,1);
        % Find closest point to plot vertical line
        [~, idx] = min(abs(plot_x - edge_x));
        if ~isempty(idx)
            plot([edge_x, edge_x], [min(field_with_diff(plot_indices)), max(params.fieldStrength(plot_indices))], ...
                 '--', 'Color', edge_colors(i,:), 'LineWidth', 1.5);
            
            % Add annotation for edge loss
            text_y = mean([min(field_with_diff(plot_indices)), max(params.fieldStrength(plot_indices))]);
            text(edge_x, text_y, sprintf('  %.1f dB', diffraction_info.edge_losses(i)), ...
                 'FontWeight', 'bold', 'Color', edge_colors(i,:));
        end
    end
    
    xlabel('Distance (m)', 'FontWeight', 'bold');
    ylabel('Field Strength (dB)', 'FontWeight', 'bold');
    title(sprintf('Field Strength Before and After Diffraction (Loss = %.2f dB, %s Method)', ...
          diffraction_info.total_loss, upper(params.method)), 'FontWeight', 'bold');
    legend('Original Field', 'With Diffraction', 'Location', 'best');
    grid on;
    
    % Add timestamp and user information
    annotation('textbox', [0.01, 0.01, 0.6, 0.03], 'String', ...
        'Current Date and Time (UTC - YYYY-MM-DD HH:MM:SS formatted): 2025-06-01 15:54:12, User: DAYALOKESH', ...
        'EdgeColor', 'none', 'FontSize', 8);
    
    % Save figure if output file specified
    if ~isempty(outputFile)
        [filepath, name, ~] = fileparts(outputFile);
        figfile = fullfile(filepath, [name, '_plot.fig']);
        saveas(gcf, figfile);
        
        % Also save as PNG with high resolution
        pngfile = fullfile(filepath, [name, '_plot.png']);
        print(pngfile, '-dpng', '-r300');
        
        fprintf('Plots saved to: %s and %s\n', figfile, pngfile);
    end
end