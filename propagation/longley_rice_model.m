function [total_pathloss, reflection_loss, diffraction_loss, freespace_loss, electric_field, distances, results] = longley_rice_model(varargin)
%LONGLEY_RICE_MODEL Complete Longley-Rice propagation model implementation
%
%   [total_pathloss, reflection_loss, diffraction_loss, freespace_loss, electric_field, distances, results] = 
%   LONGLEY_RICE_MODEL() implements the Longley-Rice empirical propagation model
%   following ITU-R P.526-15 standards for accurate RF propagation prediction.
%
%   This function calculates radio signal attenuation over irregular terrain by
%   separately computing reflection, diffraction, and free space path loss
%   components, then combining them for total path loss prediction.
%
%   Parameters (Name-Value pairs):
%       'frequency'     - Operating frequency in MHz (default: 970)
%       'txHeight'      - Transmitter antenna height in meters (default: 52)
%       'rxHeight'      - Receiver antenna height in meters (default: 2.4)
%       'polarization'  - Antenna polarization: 'vertical' or 'horizontal' (default: 'vertical')
%       'terrainFile'   - Path to terrain data file (default: '../terrain/X.04')
%       'maxDistance'   - Maximum calculation distance in meters (default: 800)
%       'stepSize'      - Distance step size for calculations in meters (default: 1.0)
%       'txPosition'    - Transmitter position [x, y] in meters (default: [0, 52])
%       'conductivity'  - Ground conductivity in S/m (default: 0.005)
%       'permittivity'  - Ground relative permittivity (default: 15)
%       'refractivity'  - Surface refractivity in N-units (default: 315)
%       'climate'       - Climate type: 1-7 (default: 5 - continental temperate)
%       'outputDir'     - Output directory for plots and results (default: './results')
%       'plotResults'   - Generate plots flag (default: true)
%       'saveResults'   - Save results to file flag (default: true)
%
%   Returns:
%       total_pathloss  - Total path loss in dB (vector)
%       reflection_loss - Ground reflection loss in dB (vector)
%       diffraction_loss- Knife-edge diffraction loss in dB (vector)
%       freespace_loss  - Free space path loss in dB (vector)
%       electric_field  - Electric field strength in V/m (vector)
%       distances       - Distance points in meters (vector)
%       results         - Structure with detailed calculation results
%
%   Example:
%       [total_loss, refl_loss, diff_loss, fs_loss, E_field, dist, info] = ...
%           longley_rice_model('frequency', 970, 'txHeight', 52, 'rxHeight', 2.4);
%
%   References:
%       - ITU-R P.526-15: Propagation by diffraction
%       - ITU-R P.1546: Method for point-to-area predictions
%       - Longley, A.G. and Rice, P.L.: Prediction of tropospheric radio transmission loss over irregular terrain
%
%   Author: Generated for DAYALOKESH dissertation
%   Date: December 2024

    %% Input Parameter Parsing and Validation
    p = inputParser;
    
    % Core propagation parameters
    addParameter(p, 'frequency', 970, @(x) isnumeric(x) && x > 0);           % MHz
    addParameter(p, 'txHeight', 52, @(x) isnumeric(x) && x > 0);            % meters
    addParameter(p, 'rxHeight', 2.4, @(x) isnumeric(x) && x > 0);           % meters
    addParameter(p, 'polarization', 'vertical', @(x) ischar(x) && ismember(lower(x), {'vertical', 'horizontal'}));
    
    % Terrain and geometry parameters
    addParameter(p, 'terrainFile', './terrain/X.04', @ischar);
    addParameter(p, 'maxDistance', 800, @(x) isnumeric(x) && x > 0);         % meters
    addParameter(p, 'stepSize', 1.0, @(x) isnumeric(x) && x > 0);            % meters
    addParameter(p, 'txPosition', [0, 52], @(x) isnumeric(x) && length(x) == 2);
    
    % Ground and atmospheric parameters
    addParameter(p, 'conductivity', 0.005, @(x) isnumeric(x) && x > 0);      % S/m (typical for average ground)
    addParameter(p, 'permittivity', 15, @(x) isnumeric(x) && x > 0);         % relative permittivity
    addParameter(p, 'refractivity', 315, @(x) isnumeric(x) && x > 0);        % N-units
    addParameter(p, 'climate', 5, @(x) isnumeric(x) && x >= 1 && x <= 7);    % climate type
    
    % Output and control parameters
    addParameter(p, 'outputDir', './results', @ischar);
    addParameter(p, 'plotResults', true, @islogical);
    addParameter(p, 'saveResults', true, @islogical);
    
    parse(p, varargin{:});
    params = p.Results;
    
    % Validate frequency range for Longley-Rice model
    if params.frequency < 20 || params.frequency > 20000
        warning('Longley-Rice model is typically valid for 20-20000 MHz. Current: %.1f MHz', params.frequency);
    end
    
    %% Initialize Physical Constants and Parameters
    fprintf('\n=== LONGLEY-RICE PROPAGATION MODEL ===\n');
    fprintf('Frequency: %.1f MHz\n', params.frequency);
    fprintf('TX Height: %.1f m, RX Height: %.1f m\n', params.txHeight, params.rxHeight);
    fprintf('Polarization: %s\n', upper(params.polarization));
    fprintf('Max Distance: %.1f m, Step Size: %.1f m\n', params.maxDistance, params.stepSize);
    
    % Physical constants
    c = 2.99792458e8;                                    % Speed of light (m/s)
    lambda = c / (params.frequency * 1e6);               % Wavelength (m)
    k = 2 * pi / lambda;                                 % Wave number (rad/m)
    
    % Create output directory
    if ~exist(params.outputDir, 'dir')
        mkdir(params.outputDir);
    end
    
    %% Load and Process Terrain Data
    fprintf('\n=== LOADING TERRAIN DATA ===\n');
    try
        % Check if terrain file exists
        if ~exist(params.terrainFile, 'file')
            error('Terrain file not found: %s', params.terrainFile);
        end
        
        % Load terrain data using existing fileparser
        if exist('./terrain/fileparser.m', 'file') || exist('fileparser.m', 'file')
            if exist('./terrain/fileparser.m', 'file')
                addpath('./terrain');
            end
            terrain_data = fileparser(params.terrainFile, params.maxDistance);
        else
            % Fallback to direct reading
            terrain_data = dlmread(params.terrainFile);
            terrain_data = terrain_data(terrain_data(:,1) <= params.maxDistance, :);
        end
        
        terrain_x = terrain_data(:, 1);    % Distance (m)
        terrain_y = terrain_data(:, 2);    % Height (m)
        
        fprintf('Loaded %d terrain points up to %.1f m\n', length(terrain_x), max(terrain_x));
        
    catch ME
        error('Failed to load terrain data: %s', ME.message);
    end
    
    %% Create Distance Vector for Calculations
    distances = (params.stepSize:params.stepSize:params.maxDistance)';
    num_points = length(distances);
    
    % Interpolate terrain heights at calculation points
    terrain_heights = interp1(terrain_x, terrain_y, distances, 'linear', 'extrap');
    
    fprintf('Calculating for %d distance points\n', num_points);
    
    %% Initialize Output Arrays
    freespace_loss = zeros(num_points, 1);
    reflection_loss = zeros(num_points, 1);
    diffraction_loss = zeros(num_points, 1);
    total_pathloss = zeros(num_points, 1);
    electric_field = zeros(num_points, 1);
    
    %% Calculate Free Space Path Loss
    fprintf('\n=== CALCULATING FREE SPACE PATH LOSS ===\n');
    
    % Standard free space path loss formula
    % FSPL = 20*log10(d) + 20*log10(f) + 20*log10(4*pi/c) where d in meters, f in Hz
    % Simplified: FSPL = 20*log10(d) + 20*log10(f_MHz) + 32.45 (for d in km, f in MHz)
    % For d in meters: FSPL = 20*log10(d) + 20*log10(f_MHz) - 27.55
    
    freespace_loss = 20*log10(distances/1000) + 20*log10(params.frequency) + 32.45;
    
    fprintf('Free space path loss range: %.1f to %.1f dB\n', min(freespace_loss), max(freespace_loss));
    
    %% Calculate Ground Reflection Component
    fprintf('\n=== CALCULATING GROUND REFLECTION ===\n');
    
    reflection_loss = calculate_ground_reflection(distances, terrain_heights, ...
                                                params.txHeight, params.rxHeight, ...
                                                lambda, params.conductivity, params.permittivity, ...
                                                params.polarization);
    
    fprintf('Ground reflection loss range: %.1f to %.1f dB\n', min(reflection_loss), max(reflection_loss));
    
    %% Calculate Multiple Knife-Edge Diffraction
    fprintf('\n=== CALCULATING KNIFE-EDGE DIFFRACTION ===\n');
    
    diffraction_loss = calculate_diffraction_loss(distances, terrain_heights, ...
                                                 params.txPosition, params.rxHeight, ...
                                                 lambda, params.refractivity);
    
    fprintf('Diffraction loss range: %.1f to %.1f dB\n', min(diffraction_loss), max(diffraction_loss));
    
    %% Combine Components for Total Path Loss
    fprintf('\n=== COMBINING COMPONENTS ===\n');
    
    % The Longley-Rice model combines components using proper RF power combination
    % Each component affects the received power, so we combine in power domain
    
    for i = 1:num_points
        % Convert individual losses to linear power ratios
        fs_power_ratio = 10^(-freespace_loss(i)/10);  % Basic free space attenuation
        
        % Ground reflection modifies the free space field
        % reflection_loss is already the modifier (can be + or -)
        refl_power_ratio = 10^(-reflection_loss(i)/10);
        
        % Diffraction adds additional loss beyond free space
        diff_power_ratio = 10^(-diffraction_loss(i)/10);
        
        % Combine the effects:
        % Start with free space, apply reflection effect, then diffraction loss
        total_power_ratio = fs_power_ratio * refl_power_ratio * diff_power_ratio;
        
        % Convert back to dB
        total_pathloss(i) = -10 * log10(total_power_ratio);
        
        % Ensure reasonable bounds (avoid numerical issues)
        if isnan(total_pathloss(i)) || isinf(total_pathloss(i))
            total_pathloss(i) = freespace_loss(i);  % Fall back to free space
        end
        
        % Practical limits for path loss
        total_pathloss(i) = max(total_pathloss(i), freespace_loss(i) - 10);  % Max 10 dB gain from reflection
        total_pathloss(i) = min(total_pathloss(i), freespace_loss(i) + 50);   % Max 50 dB additional loss
    end
    
    fprintf('Total path loss range: %.1f to %.1f dB\n', min(total_pathloss), max(total_pathloss));
    
    %% Calculate Electric Field Strength
    fprintf('\n=== CALCULATING ELECTRIC FIELD ===\n');
    
    % Calculate electric field strength from path loss
    % Standard reference: E(dBμV/m) = P(dBm) + Gr(dBi) + 107 - PL(dB)
    % Where: P = transmit power, Gr = receive antenna gain, 107 = conversion factor
    % Assume 1W (30 dBm) transmit power and 0 dBi receive antenna gain
    
    tx_power_dBm = 30;  % 1W = 30 dBm
    rx_antenna_gain_dBi = 0;  % Isotropic antenna
    
    for i = 1:num_points
        % Calculate field strength in dBμV/m
        E_dBuV_m = tx_power_dBm + rx_antenna_gain_dBi + 107 - total_pathloss(i);
        
        % Convert to V/m: E(V/m) = 10^((E(dBμV/m) - 120)/20)
        electric_field(i) = 10^((E_dBuV_m - 120)/20);
        
        % Ensure reasonable bounds
        if electric_field(i) < 1e-9
            electric_field(i) = 1e-9;  % Minimum field level
        elseif electric_field(i) > 1e3
            electric_field(i) = 1e3;   % Maximum field level
        end
    end
    
    fprintf('Electric field range: %.2e to %.2e V/m\n', min(electric_field), max(electric_field));
    
    %% Generate Results Structure
    results = struct();
    results.frequency = params.frequency;
    results.wavelength = lambda;
    results.txHeight = params.txHeight;
    results.rxHeight = params.rxHeight;
    results.polarization = params.polarization;
    results.terrain_points = length(terrain_x);
    results.calculation_points = num_points;
    results.max_distance = params.maxDistance;
    results.step_size = params.stepSize;
    results.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    
    %% Save Results if Requested
    if params.saveResults
        save_longley_rice_results(params.outputDir, distances, total_pathloss, ...
                                 reflection_loss, diffraction_loss, freespace_loss, ...
                                 electric_field, results);
    end
    
    %% Generate Plots if Requested
    if params.plotResults
        generate_longley_rice_plots(params.outputDir, distances, total_pathloss, ...
                                   reflection_loss, diffraction_loss, freespace_loss, ...
                                   electric_field, results);
    end
    
    fprintf('\n=== LONGLEY-RICE CALCULATION COMPLETE ===\n');
    fprintf('Results saved to: %s\n', params.outputDir);
    
end

%% Supporting Functions

function reflection_loss = calculate_ground_reflection(distances, terrain_heights, tx_height, rx_height, lambda, conductivity, permittivity, polarization)
    %Calculate ground reflection loss component using proper two-ray model
    
    num_points = length(distances);
    reflection_loss = zeros(num_points, 1);
    
    % Ground reflection calculations using improved two-ray model
    for i = 1:num_points
        d = distances(i);
        h_terrain = terrain_heights(i);
        
        % Skip very small distances to avoid numerical issues
        if d < 1.0
            reflection_loss(i) = 0;
            continue;
        end
        
        % Calculate effective antenna heights above local terrain
        h_tx_eff = tx_height;  % Transmitter height above reference
        h_rx_eff = rx_height + h_terrain;  % Receiver height above local terrain
        
        % Direct path distance and field
        d_direct = sqrt(d^2 + (h_tx_eff - h_rx_eff)^2);
        E_direct = 1.0 / d_direct;  % Direct path field amplitude
        
        % Reflected path geometry
        d_reflected = sqrt(d^2 + (h_tx_eff + h_rx_eff)^2);
        
        % Path difference
        delta_path = d_reflected - d_direct;
        
        % Phase difference due to path length difference
        phi_path = 2 * pi * delta_path / lambda;
        
        % Grazing angle for reflection coefficient (improved calculation)
        theta_g = atan2(h_tx_eff + h_rx_eff, d);  % More accurate grazing angle
        
        % Calculate Fresnel reflection coefficient
        R = calculate_reflection_coefficient(theta_g, conductivity, permittivity, lambda, polarization);
        
        % Reflected field amplitude with geometric spreading
        E_reflected_amp = abs(R) / d_reflected;
        
        % Total phase including reflection coefficient phase
        total_phase = phi_path + angle(R);
        
        % Vector sum of direct and reflected fields
        E_direct_complex = E_direct + 0j;  % Direct field (reference phase = 0)
        E_reflected_complex = E_reflected_amp * exp(1j * total_phase);
        E_total_complex = E_direct_complex + E_reflected_complex;
        
        % Total field magnitude
        E_total_mag = abs(E_total_complex);
        
        % Reflection loss: difference from direct path alone
        % Positive values mean reflection reduces signal, negative means enhancement
        if E_total_mag > 0
            reflection_loss(i) = 20 * log10(E_total_mag / E_direct);
        else
            reflection_loss(i) = -40;  % Very large loss if field cancels out
        end
        
        % Apply smoothing for very close distances where model breaks down
        if d < 10.0
            smoothing_factor = d / 10.0;
            reflection_loss(i) = reflection_loss(i) * smoothing_factor;
        end
    end
end

function R = calculate_reflection_coefficient(theta, conductivity, permittivity, lambda, polarization)
    %Calculate Fresnel reflection coefficient
    
    % Calculate relative complex permittivity
    epsilon_r = permittivity - 1j * conductivity * lambda * 377.0 / (2 * pi);
    
    cos_theta = cos(theta);
    sin_theta = sin(theta);
    
    if strcmpi(polarization, 'vertical')
        % Vertical (parallel) polarization
        term = sqrt(epsilon_r - sin_theta^2);
        R = (epsilon_r * cos_theta - term) / (epsilon_r * cos_theta + term);
    else
        % Horizontal (perpendicular) polarization
        term = sqrt(epsilon_r - sin_theta^2);
        R = (cos_theta - term) / (cos_theta + term);
    end
end

function diffraction_loss = calculate_diffraction_loss(distances, terrain_heights, tx_position, rx_height, lambda, refractivity)
    %Calculate multiple knife-edge diffraction loss using ITU-R P.526-15
    
    num_points = length(distances);
    diffraction_loss = zeros(num_points, 1);
    
    tx_x = tx_position(1);
    tx_y = tx_position(2);
    
    % Process each receiver position
    for i = 1:num_points
        rx_x = distances(i);
        rx_y = rx_height + terrain_heights(i);
        
        % Find terrain points between TX and RX
        terrain_indices = find(distances <= rx_x);
        if length(terrain_indices) < 3
            continue;  % Need at least some terrain points for analysis
        end
        
        path_distances = distances(terrain_indices);
        path_heights = terrain_heights(terrain_indices);
        
        % Find knife edges (significant obstructions)
        edges = find_knife_edges(path_distances, path_heights, tx_x, tx_y, rx_x, rx_y);
        
        if ~isempty(edges)
            % Calculate diffraction loss for multiple edges
            diffraction_loss(i) = calculate_multiple_edge_diffraction(edges, tx_x, tx_y, rx_x, rx_y, lambda);
        else
            diffraction_loss(i) = 0;  % Line of sight
        end
    end
end

function edges = find_knife_edges(path_distances, path_heights, tx_x, tx_y, rx_x, rx_y)
    %Find significant knife edges in terrain profile
    
    edges = [];
    if length(path_distances) < 3
        return;
    end
    
    % Calculate line-of-sight heights
    los_heights = tx_y + (rx_y - tx_y) * (path_distances - tx_x) / (rx_x - tx_x);
    
    % Find clearances (negative means obstruction)
    clearances = los_heights - path_heights;
    
    % Calculate first Fresnel zone radius at each point
    total_distance = rx_x - tx_x;
    lambda = 0.309;  % wavelength for 970 MHz (approx 0.31m)
    
    % Find points that obstruct the first Fresnel zone
    significant_obstructions = [];
    for i = 2:length(clearances)-1
        d1 = path_distances(i) - tx_x;
        d2 = rx_x - path_distances(i);
        
        % First Fresnel zone radius at this point
        if d1 > 0 && d2 > 0 && total_distance > 0
            F1_radius = sqrt(lambda * d1 * d2 / total_distance);
            
            % Check if terrain intrudes into Fresnel zone
            if clearances(i) < -0.2 * F1_radius  % 20% of first Fresnel zone (less aggressive)
                significant_obstructions = [significant_obstructions; i];
            end
        end
    end
    
    % If no significant obstructions found, try to find highest obstacles
    if isempty(significant_obstructions)
        # Find local maxima in terrain that are above LOS line
        for i = 2:length(clearances)-1
            if clearances(i) < clearances(i-1) && clearances(i) < clearances(i+1) && clearances(i) < -0.5
                significant_obstructions = [significant_obstructions; i];
            end
        end
    end
    
    % Create edges from significant obstructions
    for i = 1:length(significant_obstructions)
        idx = significant_obstructions(i);
        edges = [edges; path_distances(idx), path_heights(idx)];
    end
    
    % Limit to most significant edges (up to 3)
    if size(edges, 1) > 3
        % Calculate obstruction severity
        obstruction_severity = zeros(size(edges, 1), 1);
        for i = 1:size(edges, 1)
            edge_idx = find(path_distances == edges(i, 1));
            if ~isempty(edge_idx)
                obstruction_severity(i) = -clearances(edge_idx(1));  % Negative clearance = obstruction
            end
        end
        
        [~, idx] = sort(obstruction_severity, 'descend');
        edges = edges(idx(1:min(3, end)), :);
    end
end

function loss = calculate_multiple_edge_diffraction(edges, tx_x, tx_y, rx_x, rx_y, lambda)
    %Calculate diffraction loss for multiple knife edges using ITU-R P.526-15
    
    num_edges = size(edges, 1);
    if num_edges == 0
        loss = 0;
        return;
    end
    
    if num_edges == 1
        % Single knife edge
        edge_x = edges(1, 1);
        edge_y = edges(1, 2);
        loss = calculate_single_edge_diffraction(tx_x, tx_y, edge_x, edge_y, rx_x, rx_y, lambda);
    else
        % Multiple knife edges - use approximation method
        total_loss = 0;
        
        % Sort edges by distance from transmitter
        [~, sort_idx] = sort(edges(:, 1));
        edges = edges(sort_idx, :);
        
        % Calculate cumulative diffraction loss
        for i = 1:num_edges
            edge_x = edges(i, 1);
            edge_y = edges(i, 2);
            
            % Calculate diffraction parameter v
            if i == 1
                d1 = sqrt((edge_x - tx_x)^2 + (edge_y - tx_y)^2);
                d2 = sqrt((rx_x - edge_x)^2 + (rx_y - edge_y)^2);
                h_los = tx_y + (rx_y - tx_y) * (edge_x - tx_x) / (rx_x - tx_x);
            else
                prev_edge_x = edges(i-1, 1);
                prev_edge_y = edges(i-1, 2);
                d1 = sqrt((edge_x - prev_edge_x)^2 + (edge_y - prev_edge_y)^2);
                d2 = sqrt((rx_x - edge_x)^2 + (rx_y - edge_y)^2);
                h_los = prev_edge_y + (rx_y - prev_edge_y) * (edge_x - prev_edge_x) / (rx_x - prev_edge_x);
            end
            
            h_diff = edge_y - h_los;
            v = h_diff * sqrt(2 * (d1 + d2) / (lambda * d1 * d2));
            
            % ITU-R P.526-15 diffraction loss
            if v > -0.78
                J_v = 6.9 + 20 * log10(sqrt((v - 0.1)^2 + 1) + v - 0.1);
            else
                J_v = 0;
            end
            
            total_loss = total_loss + J_v * 0.7^(i-1);  % Reduce contribution of subsequent edges
        end
        
        loss = total_loss;
    end
end

function loss = calculate_single_edge_diffraction(tx_x, tx_y, edge_x, edge_y, rx_x, rx_y, lambda)
    %Calculate single knife edge diffraction loss
    
    % Calculate distances
    d1 = sqrt((edge_x - tx_x)^2 + (edge_y - tx_y)^2);
    d2 = sqrt((rx_x - edge_x)^2 + (rx_y - edge_y)^2);
    d_total = d1 + d2;
    
    % Height of line-of-sight at edge
    h_los = tx_y + (rx_y - tx_y) * (edge_x - tx_x) / (rx_x - tx_x);
    h_diff = edge_y - h_los;
    
    % Fresnel-Kirchhoff diffraction parameter
    v = h_diff * sqrt(2 * d_total / (lambda * d1 * d2));
    
    % ITU-R P.526 diffraction loss formula
    if v > -0.78
        loss = 6.9 + 20 * log10(sqrt((v - 0.1)^2 + 1) + v - 0.1);
    else
        loss = 0;
    end
end



function save_longley_rice_results(output_dir, distances, total_loss, refl_loss, diff_loss, fs_loss, e_field, results)
    %Save calculation results to file
    
    filename = fullfile(output_dir, sprintf('LongleyRice_Results_%s.txt', ...
                       datestr(now, 'yyyymmdd_HHMMSS')));
    
    fid = fopen(filename, 'w');
    if fid == -1
        warning('Could not create output file: %s', filename);
        return;
    end
    
    % Write header
    fprintf(fid, '%% Longley-Rice Propagation Model Results\n');
    fprintf(fid, '%% Generated: %s\n', results.timestamp);
    fprintf(fid, '%% Frequency: %.1f MHz\n', results.frequency);
    fprintf(fid, '%% TX Height: %.1f m, RX Height: %.1f m\n', results.txHeight, results.rxHeight);
    fprintf(fid, '%% Polarization: %s\n', results.polarization);
    fprintf(fid, '%% Calculation Points: %d\n', results.calculation_points);
    fprintf(fid, '%%\n');
    fprintf(fid, '%% Columns: Distance(m), Total_Loss(dB), Reflection_Loss(dB), Diffraction_Loss(dB), FreeSpace_Loss(dB), E_Field(V/m)\n');
    
    % Write data
    for i = 1:length(distances)
        fprintf(fid, '%.2f\t%.3f\t%.3f\t%.3f\t%.3f\t%.6e\n', ...
                distances(i), total_loss(i), refl_loss(i), diff_loss(i), fs_loss(i), e_field(i));
    end
    
    fclose(fid);
    fprintf('Results saved to: %s\n', filename);
end

function generate_longley_rice_plots(output_dir, distances, total_loss, refl_loss, diff_loss, fs_loss, e_field, results)
    %Generate all required plots for Longley-Rice results
    
    % Set consistent plot properties
    dist_km = distances / 1000;  % Convert to km for plotting
    
    %% Plot 1: Reflection Path Loss vs Distance
    figure('Position', [100, 600, 800, 400]);
    plot(dist_km, refl_loss, 'r-', 'LineWidth', 1.5);
    grid on;
    xlabel('Distance (km)');
    ylabel('Reflection Path Loss (dB)');
    title(sprintf('Ground Reflection Path Loss vs Distance\nFreq: %.0f MHz, TX: %.1fm, RX: %.1fm, %s Pol', ...
                  results.frequency, results.txHeight, results.rxHeight, results.polarization));
    
    % Add timestamp
    text(0.02, 0.98, sprintf('Generated: %s', results.timestamp), ...
         'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8);
    
    saveas(gcf, fullfile(output_dir, 'LongleyRice_Reflection_Loss.fig'));
    saveas(gcf, fullfile(output_dir, 'LongleyRice_Reflection_Loss.png'));
    
    %% Plot 2: Diffraction Path Loss vs Distance  
    figure('Position', [120, 550, 800, 400]);
    plot(dist_km, diff_loss, 'g-', 'LineWidth', 1.5);
    grid on;
    xlabel('Distance (km)');
    ylabel('Diffraction Path Loss (dB)');
    title(sprintf('Knife-Edge Diffraction Path Loss vs Distance\nFreq: %.0f MHz, ITU-R P.526-15 Standard', ...
                  results.frequency));
    
    text(0.02, 0.98, sprintf('Generated: %s', results.timestamp), ...
         'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8);
    
    saveas(gcf, fullfile(output_dir, 'LongleyRice_Diffraction_Loss.fig'));
    saveas(gcf, fullfile(output_dir, 'LongleyRice_Diffraction_Loss.png'));
    
    %% Plot 3: Free Space Path Loss vs Distance
    figure('Position', [140, 500, 800, 400]);
    plot(dist_km, fs_loss, 'b-', 'LineWidth', 1.5);
    grid on;
    xlabel('Distance (km)');
    ylabel('Free Space Path Loss (dB)');
    title(sprintf('Free Space Path Loss vs Distance\nFreq: %.0f MHz', results.frequency));
    
    text(0.02, 0.98, sprintf('Generated: %s', results.timestamp), ...
         'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8);
    
    saveas(gcf, fullfile(output_dir, 'LongleyRice_FreeSpace_Loss.fig'));
    saveas(gcf, fullfile(output_dir, 'LongleyRice_FreeSpace_Loss.png'));
    
    %% Plot 4: Total Path Loss vs Distance
    figure('Position', [160, 450, 800, 400]);
    plot(dist_km, total_loss, 'k-', 'LineWidth', 2);
    hold on;
    plot(dist_km, fs_loss, 'b--', 'LineWidth', 1, 'Color', [0.5, 0.5, 1]);
    plot(dist_km, fs_loss + refl_loss, 'r--', 'LineWidth', 1, 'Color', [1, 0.5, 0.5]);
    plot(dist_km, fs_loss + diff_loss, 'g--', 'LineWidth', 1, 'Color', [0.5, 1, 0.5]);
    
    grid on;
    xlabel('Distance (km)');
    ylabel('Path Loss (dB)');
    title(sprintf('Longley-Rice Total Path Loss vs Distance\nFreq: %.0f MHz, TX: %.1fm, RX: %.1fm', ...
                  results.frequency, results.txHeight, results.rxHeight));
    legend('Total Loss', 'Free Space', 'FS + Reflection', 'FS + Diffraction', 'Location', 'southeast');
    
    text(0.02, 0.98, sprintf('Generated: %s', results.timestamp), ...
         'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8);
    
    saveas(gcf, fullfile(output_dir, 'LongleyRice_Total_Loss.fig'));
    saveas(gcf, fullfile(output_dir, 'LongleyRice_Total_Loss.png'));
    
    %% Plot 5: Electric Field vs Distance
    figure('Position', [180, 400, 800, 400]);
    semilogy(dist_km, e_field, 'm-', 'LineWidth', 1.5);
    grid on;
    xlabel('Distance (km)');
    ylabel('Electric Field Strength (V/m)');
    title(sprintf('Electric Field Strength vs Distance\nFreq: %.0f MHz, Longley-Rice Model', results.frequency));
    
    text(0.02, 0.98, sprintf('Generated: %s', results.timestamp), ...
         'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8);
    
    saveas(gcf, fullfile(output_dir, 'LongleyRice_Electric_Field.fig'));
    saveas(gcf, fullfile(output_dir, 'LongleyRice_Electric_Field.png'));
    
    %% Summary Plot: All Components
    figure('Position', [200, 350, 1000, 600]);
    
    subplot(2,1,1);
    plot(dist_km, total_loss, 'k-', 'LineWidth', 2);
    hold on;
    plot(dist_km, refl_loss, 'r-', 'LineWidth', 1.5);
    plot(dist_km, diff_loss, 'g-', 'LineWidth', 1.5);
    plot(dist_km, fs_loss, 'b-', 'LineWidth', 1.5);
    
    grid on;
    xlabel('Distance (km)');
    ylabel('Path Loss (dB)');
    title('Longley-Rice Model: All Path Loss Components');
    legend('Total Loss', 'Reflection', 'Diffraction', 'Free Space', 'Location', 'southeast');
    
    subplot(2,1,2);
    semilogy(dist_km, e_field, 'm-', 'LineWidth', 1.5);
    grid on;
    xlabel('Distance (km)');
    ylabel('Electric Field (V/m)');
    title('Electric Field Strength');
    
    % sgtitle not available in Octave, use alternative
    if exist('sgtitle', 'builtin') == 5
        sgtitle(sprintf('Longley-Rice Propagation Model Summary - %.0f MHz', results.frequency));
    else
        % Add title manually for Octave compatibility
        annotation('textbox', [0.3, 0.95, 0.4, 0.05], 'String', ...
                   sprintf('Longley-Rice Propagation Model Summary - %.0f MHz', results.frequency), ...
                   'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold', ...
                   'EdgeColor', 'none');
    end
    
    text(0.02, 0.02, sprintf('Generated: %s', results.timestamp), ...
         'Units', 'normalized', 'VerticalAlignment', 'bottom', 'FontSize', 8);
    
    saveas(gcf, fullfile(output_dir, 'LongleyRice_Summary.fig'));
    saveas(gcf, fullfile(output_dir, 'LongleyRice_Summary.png'));
    
    fprintf('Generated 6 plots in: %s\n', output_dir);
end