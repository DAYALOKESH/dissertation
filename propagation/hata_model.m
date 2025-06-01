% Add utils directory to path to access physical_constants
if ~isdeployed
    [current_dir, ~, ~] = fileparts(mfilename('fullpath'));
    utils_path = fullfile(current_dir, '..', 'utils');
    if exist(utils_path, 'dir') && ~ismember(utils_path, strsplit(path, pathsep))
        addpath(utils_path);
    end
end

% hata_model: Implements the Urban Hata propagation model to predict path loss
% and electric field strength.
%
% The Hata model is an empirical formulation based on measurements made in Tokyo, Japan.
% This function specifically implements the Urban Hata model, which is suitable for
% built-up city environments.
%
% Args:
%   varargin: Name-value pairs for model parameters (see inputParser setup).
%
% Returns:
%   path_loss (double array): Path loss in dB at specified distances.
%   field_dB (double array): Electric field strength in dBμV/m.
%   E_field (double array): Electric field strength in V/m.
function [path_loss, field_dB, E_field] = hata_model(varargin)
    phys_const = physical_constants();

    % Hata Model specific constants for Urban environment (Okumura-Hata model)
    H_C1 = 69.55; % Constant offset for path loss in dB
    H_C2 = 26.16; % Factor for log10(frequency)
    H_C3 = 13.82; % Factor for log10(txHeight)
    H_C4 = 44.9;  % Constant in distance-dependent term coefficient
    H_C5 = 6.55;  % Factor for log10(txHeight) in distance-dependent term coefficient

    % Constants for mobile station antenna height correction factor (a_hm)
    % This correction is for large cities.
    AHM_C1 = 3.2;   % Coefficient for the squared log term
    AHM_C2 = 11.75; % Multiplier for rxHeight inside log10
    AHM_C3 = 4.97;  % Constant offset for a_hm

    % Input Parser: Manages and validates input parameters for the function.
    % This allows for flexible input via name-value pairs and sets default values.
    p = inputParser;
    p.addParameter('frequency', 970, @isnumeric);         % MHz (Operating frequency)
    p.addParameter('distances', 0:10:700, @isnumeric);    % meters
    p.addParameter('txHeight', 50, @isnumeric);           % Base station height (m)
    p.addParameter('rxHeight', 1.5, @isnumeric);          % Mobile station height (m)
    p.addParameter('txPower', 43, @isnumeric);            % dBm (43 dBm = 20W)
    p.addParameter('normalizeIdx', [], @(x) isempty(x) || isnumeric(x)); % Idx for normalization
    p.addParameter('refFieldDB', 0, @isnumeric);          % Reference field value for normalization
    p.addParameter('outputFile', '', @ischar);
    p.addParameter('plotResults', false, @islogical);
    p.parse(varargin{:});
    
    params = p.Results;
    
    % --- Parameter Validation Checks for Hata Model ---
    % Frequency (params.frequency): Valid range: 150 MHz to 1500 MHz.
    if params.frequency < 150 || params.frequency > 1500
        warning('Hata model validity range for frequency is 150-1500 MHz. Current frequency: %.1f MHz might lead to inaccurate results.', params.frequency);
    end

    % Base Station Height (params.txHeight): Valid range: 30 m to 200 m.
    if params.txHeight < 30 || params.txHeight > 200
        warning('Hata model validity range for base station height (txHeight) is 30-200 m. Current height: %.1f m might lead to inaccurate results.', params.txHeight);
    end

    % Mobile Station Height (params.rxHeight): Valid range: 1 m to 10 m.
    if params.rxHeight < 1 || params.rxHeight > 10
        warning('Hata model validity range for mobile station height (rxHeight) is 1-10 m. Current height: %.1f m might lead to inaccurate results.', params.rxHeight);
    end

    % Convert distances to km (Hata model uses km)
    distances_m = params.distances;
    distances_km = distances_m / 1000;
    
    % Check minimum distance to avoid log of zero - this is for calculation stability
    distances_km(distances_km < 0.001) = 0.001;  % Minimum 1 meter to prevent log(0) issues.

    % Distance (distances_km): Valid range: 1 km to 20 km.
    if max(distances_km) > 20
        warning('Hata model is typically most accurate for distances up to 20 km. Maximum requested distance: %.1f km might lead to reduced accuracy.', max(distances_km));
    end
    % Check minimum distance, ensuring it's not overly close if not caught by the 0.001km rule.
    % This also implies that if all distances are < 1km, a warning should appear.
    % We check if there are any distances > 0.02km (20m) that are still < 1km.
    min_dist_for_warning = min(distances_km(distances_km > 0.02)); % find min distance among those > 20m
    if ~isempty(min_dist_for_warning) && min_dist_for_warning < 1
        warning('Hata model is typically most accurate for distances from 1 km. Minimum requested distance (>.02km): %.3f km might lead to reduced accuracy.', min_dist_for_warning);
    elseif isempty(min_dist_for_warning) && min(distances_km) < 1 % All distances are <= 0.02km
         warning('Hata model is typically most accurate for distances from 1 km. All requested distances are very short (<= 20m or adjusted to 1m). Results may be less accurate.', min(distances_km));
    end
    % --- End of Parameter Validation Checks ---

    % Pre-allocation of Output Arrays: Initialize arrays for efficiency.
    num_points = length(distances_km);
    path_loss_out_temp = zeros(1, num_points); % Temporary name to avoid conflict if needed, though direct assignment is fine.
    field_dB = zeros(1, num_points);    % Field strength in dBμV/m
    E_field = zeros(1, num_points);     % Electric field strength in V/m

    % Calculate path loss by calling the local function
    path_loss = calculate_hata_path_loss_urban(params.frequency, params.txHeight, params.rxHeight, distances_km, ...
                                               H_C1, H_C2, H_C3, H_C4, H_C5, ...
                                               AHM_C1, AHM_C2, AHM_C3);
    
    % Calculate received power in dBm from transmit power and path loss (already vectorized).
    received_power_dBm = params.txPower - path_loss;
    
    % E-field Strength Calculation (Vectorized):
    % Convert received power from dBm to electric field strength (V/m).
    c = phys_const.c;  % Speed of light (m/s) from physical_constants
    lambda = c / (params.frequency * 1e6);  % Wavelength in meters (frequency is in MHz - this is a scalar)
    
    % Convert received power from dBm to watts (vectorized):
    % received_power_dBm is a 1xN vector.
    % power_watts_vec will also be a 1xN vector.
    power_watts_vec = 10.^((received_power_dBm - 30) / 10); % Element-wise power (dBm to W)

    % Calculate E-field strength (vectorized):
    % power_watts_vec is 1xN, lambda is scalar.
    % E_field will be a 1xN vector.
    % Formula derived from P_rx = E^2 * lambda^2 / (480 * pi^2)
    % So, E = sqrt(P_rx * 480 * pi^2 / lambda^2).
    E_field = sqrt(power_watts_vec * (480 * phys_const.PI^2 / (lambda^2))); % Scalar part multiplied by vector
    
    % Convert E-field (V/m) to dB relative to 1 microVolt/m (dBμV/m) (vectorized):
    % E_field is a 1xN vector.
    % field_dB will be a 1xN vector.
    field_dB = 20 * log10(E_field * 1e6); % Multiplication by 1e6 is element-wise, log10 is element-wise.

    % Normalization of Field Strength:
    % Adjusts the calculated field strength values if a reference point and value are provided.
    % This can be used to calibrate the model output to a known measurement.
    if ~isempty(params.normalizeIdx) && params.normalizeIdx > 0 && params.normalizeIdx <= num_points
        ref_idx = params.normalizeIdx; % Index of the reference distance/point
        
        % If a specific reference field value (params.refFieldDB) is provided AND it's not zero,
        % normalize field_dB so that field_dB(ref_idx) equals params.refFieldDB.
        if params.refFieldDB ~= 0
            fprintf('Normalization: field_dB at index %d (distance %.1fm) will be set to %.2f dB.\n', ...
                    ref_idx, distances_m(ref_idx), params.refFieldDB);
            field_offset = params.refFieldDB - field_dB(ref_idx);
            field_dB = field_dB + field_offset; % Adjust all field_dB values by this offset
        else
            % If params.refFieldDB is 0 (default), this implies no specific target value was given for normalization.
            % The original 'else' branch implemented a normalization:
            %   ref_field = field_dB(ref_idx);
            %   field_dB = ref_field + (path_loss(ref_idx) - path_loss);
            % This logic makes field_dB(i) + path_loss(i) = constant for all i, meaning all points
            % would effectively have the same received power as the reference point. This is generally
            % not desired for field strength representation, as it makes the field_dB plot merely an
            % inverted and shifted path_loss plot relative to the reference point.
            % Therefore, this alternative normalization is removed.
            % If a user wishes to normalize such that the calculated field_dB(ref_idx) becomes, for example, 0 dB
            % for relative comparison, they should:
            % 1. Run the model once to observe the calculated field_dB(ref_idx).
            % 2. Re-run the model, setting params.refFieldDB = -observed_field_dB_at_ref_idx.
            fprintf('Normalization: normalizeIdx provided, but refFieldDB is 0. No field strength normalization applied based on path loss differences.\n');
            % Original problematic line was: field_dB = field_dB(ref_idx) + (path_loss(ref_idx) - path_loss); % This line is now removed.
        end
    end
    
    % Display results for a sample point (halfway through the array)
    sample_idx = min(floor(num_points/2), num_points);
    fprintf('\nUrban Hata Model Results:\n');
    fprintf('- Frequency: %.1f MHz\n', params.frequency);
    fprintf('- Base station height: %.1f m\n', params.txHeight);
    fprintf('- Mobile station height: %.1f m\n', params.rxHeight);
    fprintf('- Sample distance: %.1f m (%.3f km)\n', distances_m(sample_idx), distances_km(sample_idx));
    fprintf('- Sample path loss: %.2f dB\n', path_loss(sample_idx));
    fprintf('- Sample field strength: %.2f dB\n', field_dB(sample_idx));
    fprintf('- Sample E-field: %.6e V/m\n', E_field(sample_idx));
    
    % Save results to file if requested
    if ~isempty(params.outputFile)
        saveResults(params, distances_m, path_loss, field_dB, E_field);
    end
    
    % Plot results if requested
    if params.plotResults
        plotHataResults(params, distances_m, path_loss, field_dB, E_field);
    end
end

function saveResults(params, distances, path_loss, field_dB, E_field)
    % Save results to a file
    fid = fopen(params.outputFile, 'w');
    
    if fid == -1
        warning('Could not open output file for writing.');
        return;
    end
    
    % Write header with timestamp
    fprintf(fid, '# Urban Hata Model Propagation Results\n');
    timestamp_str = datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd HH:mm:ss');
    fprintf(fid, '# Current Date and Time (UTC - %s formatted): %s\n', 'yyyy-MM-dd HH:mm:ss', timestamp_str);
    % Write parameters
    fprintf(fid, 'Parameters:\n');
    fprintf(fid, '- Frequency: %.2f MHz\n', params.frequency);
    fprintf(fid, '- Base station height: %.2f m\n', params.txHeight);
    fprintf(fid, '- Mobile station height: %.2f m\n', params.rxHeight);
    fprintf(fid, '- Transmitter power: %.2f dBm\n\n', params.txPower);
    
    % Write results in columns
    fprintf(fid, 'Distance(m)\tPathLoss(dB)\tFieldStrength(dB)\tE-field(V/m)\n');
    for i = 1:length(distances)
        fprintf(fid, '%.6f\t%.6f\t%.6f\t%.6e\n', ...
            distances(i), path_loss(i), field_dB(i), E_field(i));
    end
    
    % Close file
    fclose(fid);
    fprintf('Results saved to: %s\n', params.outputFile);
    
    % Also save as .dat file with just distance and field strength.
    % This format might be used for compatibility with other tools or for quick plotting (e.g., EFIE comparison).
    [filepath, name, ~] = fileparts(params.outputFile);
    dat_file = fullfile(filepath, [name, '.dat']);
    fid = fopen(dat_file, 'w');
    
    if fid ~= -1
        for i = 1:length(distances)
            fprintf(fid, '%.6f\t%.6f\n', distances(i), field_dB(i));
        end
        fclose(fid);
        fprintf('Simple format saved to: %s\n', dat_file);
    end
end

function plotHataResults(params, distances, path_loss, field_dB, E_field)
    % Create figure with multiple subplots
    figure('Name', 'Urban Hata Model Results', 'Position', [100, 100, 900, 700]);
    
    % Plot 1: Path Loss
    subplot(3, 1, 1);
    plot(distances, path_loss, 'r-', 'LineWidth', 2);
    xlabel('Distance (m)', 'FontWeight', 'bold');
    ylabel('Path Loss (dB)', 'FontWeight', 'bold');
    title(['Urban Hata Model: Path Loss (f=', num2str(params.frequency), ' MHz, hb=', ...
           num2str(params.txHeight), 'm, hm=', num2str(params.rxHeight), 'm)'], 'FontWeight', 'bold');
    grid on;
    
    % Plot 2: Field Strength in dB
    subplot(3, 1, 2);
    plot(distances, field_dB, 'b-', 'LineWidth', 2);
    xlabel('Distance (m)', 'FontWeight', 'bold');
    ylabel('Field Strength (dB)', 'FontWeight', 'bold');
    title('Urban Hata Model: Field Strength', 'FontWeight', 'bold');
    grid on;
    
    % Plot 3: Electric Field Strength
    subplot(3, 1, 3);
    plot(distances, E_field, 'm-', 'LineWidth', 2);
    xlabel('Distance (m)', 'FontWeight', 'bold');
    ylabel('Electric Field (V/m)', 'FontWeight', 'bold');
    title('Urban Hata Model: Electric Field Strength', 'FontWeight', 'bold');
    grid on;
    set(gca, 'YScale', 'log');  % Use log scale for E-field
    
    % Add timestamp and user information
    timestamp_str = datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd HH:mm:ss');
    annotation_str = sprintf('Current Date and Time (UTC): %s', timestamp_str);
    annotation('textbox', [0.01, 0.01, 0.6, 0.03], ...
               'String', annotation_str, ...
               'EdgeColor', 'none', ...
               'FontSize', 8, ...
               'FitBoxToText', 'on');
    
    % Save figure if output file specified
    if ~isempty(params.outputFile)
        [filepath, name, ~] = fileparts(params.outputFile);
        figfile = fullfile(filepath, [name, '_plot.fig']);
        saveas(gcf, figfile);
        
        % Also save as PNG
        pngfile = fullfile(filepath, [name, '_plot.png']);
        saveas(gcf, pngfile);
        
        fprintf('Plots saved to: %s and %s\n', figfile, pngfile);
    end
end

% --- Local Functions ---
function path_loss_calc = calculate_hata_path_loss_urban(frequency_MHz, txHeight_m, rxHeight_m, distances_km_vec, H_C1, H_C2, H_C3, H_C4, H_C5, AHM_C1, AHM_C2, AHM_C3)
    % Calculates Hata path loss for an urban environment.
    % Inputs:
    %   frequency_MHz: Operating frequency in MHz
    %   txHeight_m: Transmitter antenna height in meters
    %   rxHeight_m: Receiver antenna height in meters
    %   distances_km_vec: Vector of distances in kilometers
    %   H_C1 to H_C5: Hata model constants for path loss formula
    %   AHM_C1 to AHM_C3: Constants for a_hm calculation

    % Mobile station antenna height correction factor (a_hm).
    % This formula is for large cities and frequencies (f) typically >= 300 MHz.
    a_hm = AHM_C1 * (log10(AHM_C2 * rxHeight_m))^2 - AHM_C3;

    % Path Loss Calculation (Vectorized):
    % Common terms (independent of distance):
    common_loss_terms = H_C1 + H_C2 * log10(frequency_MHz) - ...
                        H_C3 * log10(txHeight_m) - a_hm;

    % Distance-dependent term coefficient:
    dist_term_coeff = H_C4 - H_C5 * log10(txHeight_m);

    % Vectorized path loss calculation:
    % log10 applies element-wise to distances_km_vec.
    path_loss_calc = common_loss_terms + dist_term_coeff * log10(distances_km_vec);
end % end of calculate_hata_path_loss_urban