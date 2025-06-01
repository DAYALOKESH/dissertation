%% Initialize workspace and parameters
clear; clc; close all;
fprintf('=== EFIE/MoM and Hata Model Analysis ===\n');
fprintf('Started: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('User: DAYALOKESH\n\n');

%% Define physical constants
constants = physical_constants(); % Call physical constants function
f_MHz = constants.f / 1e6;        % Frequency in MHz for Hata model

%% Setup paths and directories
% Add required paths
addpath('./utils');
addpath('./terrain');
addpath('./efie');
addpath('./propagation');

% Create results directory if it doesn't exist
if ~exist('./results', 'dir')
    mkdir('./results');
    fprintf('Created results directory\n');
end

%% Configure MoM solver parameters
fprintf('\n=== CONFIGURING SIMULATION PARAMETERS ===\n');

config = struct();
config.terrainFile = './terrain/X.04';
config.sourceX = 0.0;
config.sourceY = 442.0;     % Source height
config.obsHeight = 2.4;     % Observation height
config.grossStep = 10.0;    % 10m step size
config.grossSteps = 70;     % 70 steps = 700m max distance
config.outputDir = './results';
config.enablePlot = true;
config.enableLog = true;
config.showProgress = true;

% Hata model parameters
hata_params = struct();
hata_params.frequency = f_MHz;
hata_params.txHeight = 50;        % Base station height (m)
hata_params.rxHeight = 1.5;       % Mobile station height (m)
hata_params.txPower = 43;         % Transmitter power (dBm)
hata_params.cityType = 'small';   % City type for model

fprintf('EFIE Parameters:\n');
fprintf('- Frequency: %.2f MHz\n', f_MHz);
fprintf('- Source position: (%.1f, %.1f) m\n', config.sourceX, config.sourceY);
fprintf('- Observation height: %.1f m\n', config.obsHeight);
fprintf('- Maximum distance: %.1f m\n', config.grossStep * config.grossSteps);

fprintf('\nHata Model Parameters:\n');
fprintf('- Base station height: %.1f m\n', hata_params.txHeight);
fprintf('- Mobile station height: %.1f m\n', hata_params.rxHeight);
fprintf('- Transmitter power: %.1f dBm\n', hata_params.txPower);
fprintf('- City type: %s\n', hata_params.cityType);

%% Load terrain data
fprintf('\n=== LOADING TERRAIN DATA ===\n');
terrain_data = fileparser(config.terrainFile, config.grossStep * config.grossSteps);
X_terrain = terrain_data(:, 1);
Y_terrain = terrain_data(:, 2);
fprintf('Loaded terrain data: %d points (up to %.1f m)\n', length(X_terrain), max(X_terrain));

%% Run MoM solver for EFIE calculations
fprintf('\n=== EXECUTING EFIE/MoM SOLVER ===\n');
tic;
[J, Et, results_mom] = momsolver(...
    'terrainFile', config.terrainFile, ...
    'sourceX', config.sourceX, ...
    'sourceY', config.sourceY, ...
    'obsHeight', config.obsHeight, ...
    'grossStep', config.grossStep, ...
    'grossSteps', config.grossSteps, ...
    'outputDir', config.outputDir, ...
    'enablePlot', false, ...  % We'll create custom plots below
    'enableLog', config.enableLog, ...
    'showProgress', config.showProgress);
mom_time = toc;

fprintf('\nEFIE/MoM solver completed in %.2f seconds\n', mom_time);
fprintf('Results timestamp: %s\n', results_mom.timestamp);
fprintf('Max field: %.4e V/m\n', results_mom.maxField);
fprintf('Max current: %.4e A/m\n', results_mom.maxCurrent);

%% Run Urban Hata Model
fprintf('\n=== EXECUTING URBAN HATA MODEL ===\n');

% Create array of distances matching EFIE calculation points
tic;
x_coords = (1:length(Et)) * constants.DeltaX;

% Run Urban Hata model
[hata_pl, hata_field_dB, hata_field_V_m] = hata_model(...
    'frequency', hata_params.frequency, ...
    'distances', x_coords, ...
    'txHeight', hata_params.txHeight, ...
    'rxHeight', hata_params.rxHeight, ...
    'txPower', hata_params.txPower, ...
    'outputFile', './results/HataModel_Urban_Results.txt', ...
    'plotResults', false);  % We'll create custom plots below
hata_time = toc;

fprintf('Urban Hata model completed in %.2f seconds\n', hata_time);
fprintf('Max field: %.2f dB\n', max(hata_field_dB));
fprintf('Min path loss: %.2f dB\n', min(hata_pl));

%% Convert EFIE results to dB for plotting
% Calculate E-field in dB relative to free space
EFIE_field_dB = zeros(size(Et));
for i = 1:length(Et)
    R_obs = sqrt((config.sourceX - x_coords(i))^2 + (config.sourceY - (Y_terrain(min(i, length(Y_terrain))) + config.obsHeight))^2);
    EFIE_field_dB(i) = 20*log10(abs(Et(i))/sqrt(R_obs) + eps);
end

%% Generate plots for EFIE and Hata results
fprintf('\n=== GENERATING VISUALIZATIONS ===\n');

% Plot the Terrain
fprintf('Plotting terrain profile...\n');
figure('Units','normalized','Position',[0.1 0.6 0.4 0.3]);
plot(X_terrain, Y_terrain, 'b-', 'LineWidth', 1.5);
hold on;
plot(config.sourceX, config.sourceY, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('Distance (m)');
ylabel('Terrain Elevation (m)');
title('Terrain Profile with Source Location');
grid on;
legend('Terrain', 'Source', 'Location', 'best');
add_timestamp_annotation();
saveas(gcf, './results/Terrain_Profile.fig');
saveas(gcf, './results/Terrain_Profile.png');

% Surface current distribution (EFIE)
fprintf('Plotting surface current distribution...\n');
figure('Units','normalized','Position',[0.1 0.3 0.4 0.4]);
plot(x_coords, abs(J), 'b-', 'LineWidth', 1.5);
xlabel('Distance Along Surface (m)');
ylabel('|J| (A/m)');
title('Surface Current Distribution (EFIE Method)');
grid on;
add_timestamp_annotation();
saveas(gcf, './results/SurfaceCurrent_EFIE.fig');
saveas(gcf, './results/SurfaceCurrent_EFIE.png');

% EFIE Electric Field
fprintf('Plotting EFIE electric field...\n');
figure('Units','normalized','Position',[0.5 0.3 0.4 0.4]);
plot(x_coords, EFIE_field_dB, 'b-', 'LineWidth', 1.5);
xlabel('Distance Along Surface (m)');
ylabel('Electric Field (dB)');
title('Electric Field Strength (EFIE Method)');
grid on;
add_timestamp_annotation();
saveas(gcf, './results/ElectricField_EFIE.fig');
saveas(gcf, './results/ElectricField_EFIE.png');

% Hata model path loss
fprintf('Plotting Hata model path loss...\n');
figure('Units','normalized','Position',[0.5 0.6 0.4 0.3]);
plot(x_coords, hata_pl, 'r-', 'LineWidth', 1.5);
xlabel('Distance (m)');
ylabel('Path Loss (dB)');
title(['Hata Urban Model Path Loss (f=', num2str(hata_params.frequency), ' MHz, hb=', ...
      num2str(hata_params.txHeight), 'm, hm=', num2str(hata_params.rxHeight), 'm)']);
grid on;
add_timestamp_annotation();
saveas(gcf, './results/HataModel_PathLoss.fig');
saveas(gcf, './results/HataModel_PathLoss.png');

% Hata model field strength
fprintf('Plotting Hata model field strength...\n');
figure('Units','normalized','Position',[0.1 0.1 0.4 0.4]);
plot(x_coords, hata_field_dB, 'r-', 'LineWidth', 1.5);
xlabel('Distance (m)');
ylabel('Field Strength (dB)');
title(['Hata Urban Model Field Strength (f=', num2str(hata_params.frequency), ' MHz)']);
grid on;
add_timestamp_annotation();
saveas(gcf, './results/HataModel_FieldStrength.fig');
saveas(gcf, './results/HataModel_FieldStrength.png');

% Hata model E-field
fprintf('Plotting Hata model E-field...\n');
figure('Units','normalized','Position',[0.5 0.1 0.4 0.4]);
plot(x_coords, hata_field_V_m, 'r-', 'LineWidth', 1.5);
xlabel('Distance (m)');
ylabel('Electric Field (V/m)');
title(['Hata Urban Model E-Field (f=', num2str(hata_params.frequency), ' MHz)']);
grid on;
set(gca, 'YScale', 'log');  % Use log scale for E-field
add_timestamp_annotation();
saveas(gcf, './results/HataModel_EField.fig');
saveas(gcf, './results/HataModel_EField.png');

%% Save summary information
fprintf('\n=== GENERATING SUMMARY REPORT ===\n');
fileID = fopen('./results/Calculation_Summary.txt', 'w');
fprintf(fileID, 'EFIE and Hata Model Calculation Summary\n');
fprintf(fileID, 'Current Date and Time (UTC): 2025-06-01 14:47:51\n');
fprintf(fileID, 'User: DAYALOKESH\n\n');

fprintf(fileID, 'Parameters:\n');
fprintf(fileID, '- Frequency: %.2f MHz\n', f_MHz);
fprintf(fileID, '- Wavelength: %.4f m\n', constants.Lambda);
fprintf(fileID, '- Source Position: (%.1f, %.1f) m\n', config.sourceX, config.sourceY);
fprintf(fileID, '- Observation Height: %.1f m\n', config.obsHeight);
fprintf(fileID, '- Base Station Height (Hata): %.1f m\n', hata_params.txHeight);
fprintf(fileID, '- Mobile Station Height (Hata): %.1f m\n\n', hata_params.rxHeight);

fprintf(fileID, 'Performance:\n');
fprintf(fileID, '- EFIE/MoM Calculation Time: %.2f seconds\n', mom_time);
fprintf(fileID, '- Hata Model Calculation Time: %.2f seconds\n', hata_time);
fprintf(fileID, '- Total Calculation Time: %.2f seconds\n\n', mom_time + hata_time);

fprintf(fileID, 'Results:\n');
fprintf(fileID, '- Max EFIE Field: %.2f dB\n', max(EFIE_field_dB));
fprintf(fileID, '- Max EFIE Current: %.4e A/m\n', max(abs(J)));
fprintf(fileID, '- Max Hata Field: %.2f dB\n', max(hata_field_dB));
fprintf(fileID, '- Min Hata Path Loss: %.2f dB\n', min(hata_pl));
fclose(fileID);

%% Cleanup and finish
fprintf('\n=== CALCULATIONS COMPLETE ===\n');
fprintf('All results saved in %s\n', config.outputDir);
fprintf('Ended: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf('Total execution time: %.2f seconds\n', mom_time + hata_time);

%% Helper function to add timestamp annotation to figures
function add_timestamp_annotation()
    annotation('textbox', [0.01, 0.01, 0.6, 0.03], 'String', ...
        'Current Date and Time (UTC): 2025-06-01 14:47:51, User: DAYALOKESH', ...
        'EdgeColor', 'none', 'FontSize', 8);
end