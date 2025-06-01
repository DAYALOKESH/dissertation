function [path_loss, field_dB, E_field] = hata_model(varargin)
    % Parse input parameters
    p = inputParser;
    p.addParameter('frequency', 970, @isnumeric);         % MHz
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
    
    % Convert distances to km (Hata model uses km)
    distances_m = params.distances;
    distances_km = distances_m / 1000;
    
    % Check minimum distance to avoid log of zero
    distances_km(distances_km < 0.001) = 0.001;  % Minimum 1 meter
    
    % Calculate mobile station antenna height correction factor
    % For large city and f >= 400 MHz
    a_hm = 3.2 * (log10(11.75 * params.rxHeight))^2 - 4.97;
    
    % Pre-allocate arrays for path loss and field strength
    num_points = length(distances_km);
    path_loss = zeros(1, num_points);
    field_dB = zeros(1, num_points);
    E_field = zeros(1, num_points);
    
    % Calculate urban path loss using Hata model
    for i = 1:num_points
        % Standard Hata urban model formula
        path_loss(i) = 69.55 + 26.16 * log10(params.frequency) - ...
                      13.82 * log10(params.txHeight) - a_hm + ...
                      (44.9 - 6.55 * log10(params.txHeight)) * log10(distances_km(i));
    end
    
    % Calculate received power
    received_power_dBm = params.txPower - path_loss;
    
    % Convert received power to electric field strength
    c = 299792458;  % Speed of light (m/s)
    lambda = c / (params.frequency * 1e6);  % Wavelength in meters
    
    for i = 1:num_points
        % Convert dBm to watts
        power_watts = 10^((received_power_dBm(i) - 30) / 10);
        
        % Convert power to field strength: E = sqrt(P * 480π²/(λ²))
        E_field(i) = sqrt(power_watts * 480 * pi^2 / (lambda^2));
        
        % Field in dB relative to reference (normally 1 μV/m)
        field_dB(i) = 20 * log10(E_field(i) * 1e6);
    end
    
    % Normalize field strength if reference point provided
    if ~isempty(params.normalizeIdx) && params.normalizeIdx > 0 && params.normalizeIdx <= num_points
        ref_idx = params.normalizeIdx;
        ref_pl = path_loss(ref_idx);
        
        % If reference field value provided, use it; otherwise use calculated value
        if params.refFieldDB ~= 0
            field_offset = params.refFieldDB - field_dB(ref_idx);
            field_dB = field_dB + field_offset;
        else
            % Alternative normalization as in the reference code
            ref_field = field_dB(ref_idx);
            field_dB = ref_field + (ref_pl - path_loss);
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
    fprintf(fid, '# Current Date and Time (UTC - YYYY-MM-DD HH:MM:SS formatted): 2025-06-01 14:33:38\n');
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
    
    % Also save as .dat file with just distance and field strength for EFIE comparison
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
    annotation('textbox', [0.01, 0.01, 0.6, 0.03], ...
           'String', 'Current Date and Time (UTC): 2025-06-01 14:33:38', ...
           'EdgeColor', 'none', ...
           'FontSize', 8, ...
           'FitBoxToText', 'on');  % Optional: fits the box tightly around the text

    
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