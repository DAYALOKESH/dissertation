%% Physics Validation Script for Longley-Rice Model
% This script validates that the improved calculations follow proper physics

clear; clc; close all;

fprintf('=== LONGLEY-RICE MODEL PHYSICS VALIDATION ===\n\n');

% Add required paths
addpath('./propagation');
addpath('./terrain');

%% Run the model with a small subset for analysis
[total_loss, refl_loss, diff_loss, fs_loss, E_field, distances, results] = ...
    longley_rice_model('frequency', 970, 'txHeight', 52, 'rxHeight', 2.4, ...
                       'maxDistance', 100, 'stepSize', 10, 'plotResults', false, 'saveResults', false);

%% Physics Validation Tests

fprintf('1. FREE SPACE PATH LOSS VALIDATION:\n');
% Free space loss should increase monotonically with distance
fs_differences = diff(fs_loss);
if all(fs_differences >= 0)
    fprintf('   ✓ Free space path loss increases monotonically\n');
else
    fprintf('   ✗ Free space path loss has non-monotonic behavior\n');
end

% Check if FSPL formula is correct: FSPL = 20*log10(d_km) + 20*log10(f_MHz) + 32.45
d_km = distances(end) / 1000;
f_MHz = 970;
expected_fspl = 20*log10(d_km) + 20*log10(f_MHz) + 32.45;
actual_fspl = fs_loss(end);
fspl_error = abs(expected_fspl - actual_fspl);
if fspl_error < 0.1
    fprintf('   ✓ Free space path loss formula correct (error: %.3f dB)\n', fspl_error);
else
    fprintf('   ⚠ Free space path loss formula may have issues (error: %.3f dB)\n', fspl_error);
end

fprintf('\n2. GROUND REFLECTION VALIDATION:\n');
% Reflection loss should be reasonable (typically -10 to +10 dB for most scenarios)
max_refl = max(abs(refl_loss));
if max_refl < 15
    fprintf('   ✓ Ground reflection values are reasonable (max: %.1f dB)\n', max_refl);
else
    fprintf('   ⚠ Ground reflection values may be too large (max: %.1f dB)\n', max_refl);
end

% Check for NaN or infinite values
if any(isnan(refl_loss)) || any(isinf(refl_loss))
    fprintf('   ✗ Ground reflection contains NaN or infinite values\n');
else
    fprintf('   ✓ Ground reflection values are finite\n');
end

fprintf('\n3. DIFFRACTION VALIDATION:\n');
% Diffraction loss should be non-negative (it can only add loss, not gain)
if all(diff_loss >= 0)
    fprintf('   ✓ Diffraction loss values are non-negative\n');
else
    fprintf('   ✗ Diffraction loss has negative values\n');
end

% Check for reasonable diffraction values
max_diff = max(diff_loss);
if max_diff < 100
    fprintf('   ✓ Diffraction loss values are reasonable (max: %.1f dB)\n', max_diff);
else
    fprintf('   ⚠ Diffraction loss values may be too large (max: %.1f dB)\n', max_diff);
end

fprintf('\n4. TOTAL PATH LOSS VALIDATION:\n');
% Total loss should generally increase with distance
total_trend = polyfit(distances, total_loss, 1);  % Linear fit
if total_trend(1) > 0  % Positive slope
    fprintf('   ✓ Total path loss generally increases with distance\n');
else
    fprintf('   ⚠ Total path loss trend may be incorrect\n');
end

% Total loss should be at least as much as free space loss (minus small reflection gains)
min_reasonable_loss = fs_loss - 10;  % Allow up to 10 dB gain from reflection
if all(total_loss >= min_reasonable_loss)
    fprintf('   ✓ Total path loss is physically reasonable compared to free space\n');
else
    fprintf('   ⚠ Total path loss may have unrealistic gains\n');
end

fprintf('\n5. ELECTRIC FIELD VALIDATION:\n');
% Electric field should generally decrease with distance
E_trend = polyfit(log10(distances), log10(E_field), 1);  % Log-log fit
if E_trend(1) < 0  % Negative slope in log-log space
    fprintf('   ✓ Electric field decreases with distance\n');
else
    fprintf('   ⚠ Electric field trend may be incorrect\n');
end

% Check for reasonable electric field magnitudes
min_E = min(E_field);
max_E = max(E_field);
if min_E > 1e-12 && max_E < 1e3
    fprintf('   ✓ Electric field magnitudes are reasonable (%.2e to %.2e V/m)\n', min_E, max_E);
else
    fprintf('   ⚠ Electric field magnitudes may be unrealistic\n');
end

fprintf('\n6. CONSISTENCY VALIDATION:\n');
% Check if path loss and electric field are consistent
% E ∝ 1/sqrt(Path_Loss_Linear), so log(E) should be inversely related to Path_Loss_dB
correlation_coeff = corrcoef(total_loss, log10(E_field));
if abs(correlation_coeff(1,2)) > 0.9
    fprintf('   ✓ Path loss and electric field are strongly correlated (r = %.3f)\n', correlation_coeff(1,2));
else
    fprintf('   ⚠ Path loss and electric field correlation may be weak (r = %.3f)\n', correlation_coeff(1,2));
end

fprintf('\n=== VALIDATION SUMMARY ===\n');
fprintf('Longley-Rice model physics validation completed.\n');
fprintf('Review any warnings (⚠) or errors (✗) above.\n');
fprintf('✓ indicates passed validation tests.\n\n');

%% Display sample results for manual inspection
fprintf('SAMPLE CALCULATION RESULTS:\n');
fprintf('Distance(m)  FSPL(dB)  Refl(dB)  Diff(dB)  Total(dB)  E-field(V/m)\n');
fprintf('---------------------------------------------------------------\n');
for i = 1:min(5, length(distances))
    fprintf('%8.0f  %7.1f  %8.1f  %8.1f  %8.1f  %11.2e\n', ...
            distances(i), fs_loss(i), refl_loss(i), diff_loss(i), total_loss(i), E_field(i));
end
fprintf('   ...       ...       ...       ...       ...           ...\n');
for i = max(1, length(distances)-2):length(distances)
    fprintf('%8.0f  %7.1f  %8.1f  %8.1f  %8.1f  %11.2e\n', ...
            distances(i), fs_loss(i), refl_loss(i), diff_loss(i), total_loss(i), E_field(i));
end