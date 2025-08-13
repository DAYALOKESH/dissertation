%% Test Script for Longley-Rice Model
% This script tests the Longley-Rice propagation model implementation
% with the specified parameters from the problem statement

clear; clc; close all;
fprintf('=== LONGLEY-RICE MODEL TEST ===\n');
fprintf('Started: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

%% Add required paths
addpath('./propagation');
addpath('./terrain');
addpath('./utils');

%% Create results directory
if ~exist('./results', 'dir')
    mkdir('./results');
end

%% Set parameters as specified in problem statement
test_params = struct();
test_params.frequency = 970;        % MHz (as specified)
test_params.txHeight = 52;          % m (as specified) 
test_params.rxHeight = 2.4;         % m (as specified)
test_params.polarization = 'vertical'; % as specified
test_params.terrainFile = './terrain/X.04'; % as specified
test_params.maxDistance = 700;      % m (reasonable test distance)
test_params.stepSize = 5;           % m (for faster computation in test)
test_params.outputDir = './results';
test_params.plotResults = true;
test_params.saveResults = true;

fprintf('\nTest Parameters:\n');
fprintf('- Frequency: %d MHz\n', test_params.frequency);
fprintf('- TX Height: %.1f m\n', test_params.txHeight);
fprintf('- RX Height: %.1f m\n', test_params.rxHeight);
fprintf('- Polarization: %s\n', test_params.polarization);
fprintf('- Max Distance: %d m\n', test_params.maxDistance);
fprintf('- Step Size: %d m\n', test_params.stepSize);

%% Test the Longley-Rice model
try
    fprintf('\n=== CALLING LONGLEY-RICE MODEL ===\n');
    
    [total_pathloss, reflection_loss, diffraction_loss, freespace_loss, ...
     electric_field, distances, results] = longley_rice_model(...
        'frequency', test_params.frequency, ...
        'txHeight', test_params.txHeight, ...
        'rxHeight', test_params.rxHeight, ...
        'polarization', test_params.polarization, ...
        'terrainFile', test_params.terrainFile, ...
        'maxDistance', test_params.maxDistance, ...
        'stepSize', test_params.stepSize, ...
        'outputDir', test_params.outputDir, ...
        'plotResults', test_params.plotResults, ...
        'saveResults', test_params.saveResults);
    
    %% Display summary results
    fprintf('\n=== TEST RESULTS SUMMARY ===\n');
    fprintf('Calculation completed successfully!\n');
    fprintf('Number of distance points: %d\n', length(distances));
    fprintf('Distance range: %.1f to %.1f m\n', min(distances), max(distances));
    fprintf('\nPath Loss Results:\n');
    fprintf('- Free Space Loss: %.1f to %.1f dB\n', min(freespace_loss), max(freespace_loss));
    fprintf('- Reflection Loss: %.1f to %.1f dB\n', min(reflection_loss), max(reflection_loss));
    fprintf('- Diffraction Loss: %.1f to %.1f dB\n', min(diffraction_loss), max(diffraction_loss));
    fprintf('- Total Path Loss: %.1f to %.1f dB\n', min(total_pathloss), max(total_pathloss));
    fprintf('\nElectric Field:\n');
    fprintf('- Range: %.2e to %.2e V/m\n', min(electric_field), max(electric_field));
    
    %% Validate results
    fprintf('\n=== VALIDATION CHECKS ===\n');
    
    % Check that all outputs have same length
    if length(total_pathloss) == length(reflection_loss) && ...
       length(reflection_loss) == length(diffraction_loss) && ...
       length(diffraction_loss) == length(freespace_loss) && ...
       length(freespace_loss) == length(electric_field) && ...
       length(electric_field) == length(distances)
        fprintf('✓ All output arrays have consistent length\n');
    else
        fprintf('✗ Output arrays have inconsistent lengths\n');
    end
    
    % Check for reasonable values
    if all(freespace_loss > 0) && all(isfinite(freespace_loss))
        fprintf('✓ Free space loss values are reasonable\n');
    else
        fprintf('✗ Free space loss has invalid values\n');
    end
    
    if all(isfinite(total_pathloss)) && all(total_pathloss > 0)
        fprintf('✓ Total path loss values are reasonable\n');
    else
        fprintf('✗ Total path loss has invalid values\n');
    end
    
    if all(electric_field > 0) && all(isfinite(electric_field))
        fprintf('✓ Electric field values are reasonable\n');
    else
        fprintf('✗ Electric field has invalid values\n');
    end
    
    % Check that files were created
    expected_files = {
        'LongleyRice_Reflection_Loss.png',
        'LongleyRice_Diffraction_Loss.png', 
        'LongleyRice_FreeSpace_Loss.png',
        'LongleyRice_Total_Loss.png',
        'LongleyRice_Electric_Field.png',
        'LongleyRice_Summary.png'
    };
    
    files_created = 0;
    for i = 1:length(expected_files)
        if exist(fullfile(test_params.outputDir, expected_files{i}), 'file')
            files_created = files_created + 1;
        end
    end
    
    fprintf('✓ Created %d out of %d expected plot files\n', files_created, length(expected_files));
    
    fprintf('\n=== TEST COMPLETED SUCCESSFULLY ===\n');
    
catch ME
    fprintf('\n=== TEST FAILED ===\n');
    fprintf('Error: %s\n', ME.message);
    fprintf('Stack trace:\n');
    for i = 1:length(ME.stack)
        fprintf('  %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
    end
end

fprintf('\nTest completed at: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));