%% MATLAB Implementation for EFIE-based MoM Solver with Optimized Timing and Progress Bars
function momsolver()
    clear all; close all; clc;
    
    % ==============================
    % Initialize Timing and Logging
    % ==============================
    start_time_total = tic;
    run_timestamp = datestr(now, 'dd-mm_HH-MM-SS');
    
    % Add paths for utility functions
    addpath('../utils');
    addpath('../terrain');
    
    % Create results directory if it doesn't exist
    if ~exist('../results', 'dir')
        mkdir('../results');
    end
    
    % Initialize log file
    log_filename = sprintf('../results/run-time_%s.log', run_timestamp);
    log_file = fopen(log_filename, 'w');
    
    % Write log header
    fprintf(log_file, '====================================================\n');
    fprintf(log_file, 'MoM SOLVER EXECUTION LOG\n');
    fprintf(log_file, '====================================================\n');
    fprintf(log_file, 'Start Time: %s\n', datestr(now, 'dd/mm/yyyy HH:MM:SS'));
    fprintf(log_file, '====================================================\n\n');
    
    fprintf('Starting MoM solver execution...\n');
    
    % ==============================
    % Get Physical Constants and Setup
    % ==============================
    constants = physical_constants();
    PI = constants.PI;
    Epsilon_0 = constants.Epsilon_0;
    Mu_0 = constants.Mu_0;
    c = constants.c;
    f = constants.f;
    Lambda = constants.Lambda;
    DeltaX = constants.DeltaX;
    Omega = constants.Omega;
    Beta_0 = constants.Beta_0;
    
    % Terrain and Source Parameters
    GrossStep = 10.0;
    GrossNoSteps = 70;
    NoLinesubs = floor((GrossStep * GrossNoSteps) / DeltaX);
    
    % Source and Observation Parameters
    Xsource = 0.0;
    Ysource = 442.0;
    H = 2.4;
    
    % Read Terrain Data
    terrain_data = fileparser('../terrain/X.04');
    X_terrain = terrain_data(:, 1);
    Y_terrain = terrain_data(:, 2);
    
    % Log problem parameters
    fprintf(log_file, 'PROBLEM PARAMETERS:\n');
    fprintf(log_file, '- Matrix size: %d x %d\n', NoLinesubs, NoLinesubs);
    fprintf(log_file, '- Frequency: %.2e Hz\n', f);
    fprintf(log_file, '- Wavelength: %.4f m\n', Lambda);
    fprintf(log_file, '- Terrain points: %d\n', length(X_terrain));
    fprintf(log_file, '- Source: (%.1f, %.1f) m\n\n', Xsource, Ysource);
    
    fprintf('Matrix size: %d x %d (%.1f MB)\n', NoLinesubs, NoLinesubs, ...
        (NoLinesubs^2 * 16) / (1024^2)); % 16 bytes per complex double
    
    % ==============================
    % Build Z Matrix with Progress Bar
    % ==============================
    fprintf('\n=== Z MATRIX CONSTRUCTION ===\n');
    fprintf(log_file, 'Z MATRIX CONSTRUCTION:\n');
    
    z_matrix_start = tic;
    Z = complex(zeros(NoLinesubs, NoLinesubs));
    
    % Initialize progress bar
    fprintf('Building Z matrix...\n');
    fprintf('[');
    progress_width = 50;
    last_progress = 0;
    
    for p = 1:NoLinesubs
        % Calculate current progress
        current_progress = floor((p / NoLinesubs) * progress_width);
        
        % Update progress bar
        if current_progress > last_progress
            for i = last_progress+1:current_progress
                fprintf('=');
            end
            last_progress = current_progress;
        end
        
        % Build matrix row
        for q = 1:NoLinesubs
            if p == q
                % Self-impedance term
                R_self = segment_length(p, X_terrain, Y_terrain, DeltaX, GrossStep);
                Z(p, q) = (Beta_0^2 / (4*Omega*Epsilon_0)) * ...
                         (R_self - 1i*(2*R_self/PI)*log(1.781*Beta_0*R_self/(4*exp(1))));
            else
                % Mutual impedance term
                Rpq = distance_p_q(p, q, X_terrain, Y_terrain, DeltaX, GrossStep);
                Z(p, q) = (Beta_0^2 / (4*Omega*Epsilon_0)) * ...
                         segment_length(q, X_terrain, Y_terrain, DeltaX, GrossStep) * ...
                         (besselj(0, Beta_0*Rpq) - 1i*bessely(0, Beta_0*Rpq));
            end
        end
        
        % Display percentage every 10%
        if mod(p, floor(NoLinesubs/10)) == 0
            fprintf(' %d%%', round((p/NoLinesubs)*100));
        end
    end
    
    % Complete progress bar
    for i = last_progress+1:progress_width
        fprintf('=');
    end
    fprintf('] 100%%\n');
    
    z_matrix_time = toc(z_matrix_start);
    
    % Analyze Z matrix
    cond_Z = cond(Z);
    
    fprintf('Z matrix completed in %.2f seconds\n', z_matrix_time);
    fprintf('Condition number: %.2e', cond_Z);
    
    if cond_Z > 1e12
        fprintf(' [WARNING: ILL-CONDITIONED]\n');
        fprintf(log_file, '- Status: ILL-CONDITIONED (cond = %.2e)\n', cond_Z);
    else
        fprintf(' [OK]\n');
        fprintf(log_file, '- Status: WELL-CONDITIONED\n');
    end
    
    fprintf(log_file, '- Construction time: %.4f seconds\n', z_matrix_time);
    fprintf(log_file, '- Condition number: %.6e\n', cond_Z);
    fprintf(log_file, '- Matrix norm: %.6e\n', norm(Z, 'fro'));
    fprintf(log_file, '- Memory usage: %.1f MB\n\n', (NoLinesubs^2 * 16) / (1024^2));
    
    % ==============================
    % Build E Vector (Quick)
    % ==============================
    fprintf('\nBuilding incident field vector...');
    E = complex(zeros(NoLinesubs, 1));
    
    for p = 1:NoLinesubs
        R_source = distance_source_p(p, Xsource, Ysource, X_terrain, Y_terrain, DeltaX, GrossStep);
        E(p) = -(Beta_0^2/(4*Omega*Epsilon_0)) * ...
               (besselj(0, Beta_0*R_source) - 1i*bessely(0, Beta_0*R_source));
    end
    fprintf(' Done\n');
    
    % ==============================
    % Solve Matrix Equation with Progress Monitoring
    % ==============================
    fprintf('\n=== MATRIX SOLUTION ===\n');
    fprintf(log_file, 'MATRIX SOLUTION:\n');
    
    fprintf('Solving J = Z\\E using MATLAB backslash operator...\n');
    solve_start = tic;
    
    % Monitor memory usage before solve
    mem_info = memory;
    fprintf('Available memory: %.1f GB\n', mem_info.MemAvailableAllArrays / (1024^3));
    
    % Solve the system
    fprintf('Computing LU decomposition and back-substitution...');
    J = Z \ E;
    solve_time = toc(solve_start);
    fprintf(' Done\n');
    
    % Verify solution quality
    residual = norm(Z*J - E) / norm(E);
    fprintf('Solution completed in %.2f seconds\n', solve_time);
    fprintf('Relative residual: %.2e', residual);
    
    if residual < 1e-10
        fprintf(' [EXCELLENT]\n');
        fprintf(log_file, '- Status: EXCELLENT (residual = %.2e)\n', residual);
    elseif residual < 1e-6
        fprintf(' [GOOD]\n');
        fprintf(log_file, '- Status: GOOD (residual = %.2e)\n', residual);
    else
        fprintf(' [POOR - CHECK RESULTS]\n');
        fprintf(log_file, '- Status: POOR (residual = %.2e)\n', residual);
    end
    
    fprintf(log_file, '- Solution time: %.4f seconds\n', solve_time);
    fprintf(log_file, '- Relative residual: %.6e\n', residual);
    fprintf(log_file, '- Max current: %.6e A/m\n', max(abs(J)));
    fprintf(log_file, '- Current norm: %.6e\n\n', norm(J));
    
    % ==============================
    % Calculate Electric Field with Progress Bar
    % ==============================
    fprintf('\n=== FIELD CALCULATION ===\n');
    fprintf('Computing scattered field...\n');
    fprintf('[');
    
    Et = complex(zeros(1, NoLinesubs));
    last_progress = 0;
    
    for idx = 1:NoLinesubs
        % Update progress bar
        current_progress = floor((idx / NoLinesubs) * progress_width);
        if current_progress > last_progress
            for i = last_progress+1:current_progress
                fprintf('=');
            end
            last_progress = current_progress;
        end
        
        % Calculate scattered field
        Et_scattered = complex(0, 0);
        for n = 1:NoLinesubs
            Rn = distance_surf_obs(n, idx, X_terrain, Y_terrain, DeltaX, GrossStep, H);
            Z_val = (Beta_0^2 / (4*Omega*Epsilon_0)) * ...
                   (besselj(0, Beta_0*Rn) - 1i*bessely(0, Beta_0*Rn));
            Et_scattered = Et_scattered + J(n) * segment_length(n, X_terrain, Y_terrain, DeltaX, GrossStep) * Z_val;
        end
        
        % Add incident field
        R_obs = distance_source_obs(idx, Xsource, Ysource, X_terrain, Y_terrain, DeltaX, GrossStep, H);
        Ei = -(Beta_0^2/(4*Omega*Epsilon_0)) * ...
            (besselj(0, Beta_0*R_obs) - 1i*bessely(0, Beta_0*R_obs));
        Et(idx) = Ei - Et_scattered;
        
        % Display percentage
        if mod(idx, floor(NoLinesubs/10)) == 0
            fprintf(' %d%%', round((idx/NoLinesubs)*100));
        end
    end
    
    % Complete progress bar
    for i = last_progress+1:progress_width
        fprintf('=');
    end
    fprintf('] 100%%\n');
    
    % ==============================
    % Generate Output Files and Plots
    % ==============================
    fprintf('\n=== OUTPUT GENERATION ===\n');
    fprintf('Writing results to files...');
    write_results_timestamped(J, Et, X_terrain, Y_terrain, DeltaX, GrossStep, NoLinesubs, ...
                             Xsource, Ysource, H, run_timestamp);
    fprintf(' Done\n');
    
    fprintf('Generating plots...');
    plot_results_timestamped(J, Et, X_terrain, Y_terrain, DeltaX, GrossStep, NoLinesubs, ...
                            Xsource, Ysource, H, run_timestamp);
    fprintf(' Done\n');
    
    % ==============================
    % Final Summary
    % ==============================
    total_time = toc(start_time_total);
    
    fprintf('\n====================================================\n');
    fprintf('EXECUTION SUMMARY\n');
    fprintf('====================================================\n');
    fprintf('Z matrix construction: %8.2f s (%5.1f%%)\n', z_matrix_time, (z_matrix_time/total_time)*100);
    fprintf('Matrix solution:       %8.2f s (%5.1f%%)\n', solve_time, (solve_time/total_time)*100);
    fprintf('Other operations:      %8.2f s (%5.1f%%)\n', ...
        total_time-z_matrix_time-solve_time, ((total_time-z_matrix_time-solve_time)/total_time)*100);
    fprintf('----------------------------------------------------\n');
    fprintf('TOTAL EXECUTION TIME:  %8.2f s\n', total_time);
    fprintf('====================================================\n');
    
    fprintf(log_file, 'EXECUTION SUMMARY:\n');
    fprintf(log_file, '- Z matrix construction: %.4f seconds (%.1f%%)\n', z_matrix_time, (z_matrix_time/total_time)*100);
    fprintf(log_file, '- Matrix solution: %.4f seconds (%.1f%%)\n', solve_time, (solve_time/total_time)*100);
    fprintf(log_file, '- Total time: %.4f seconds\n', total_time);
    fprintf(log_file, '- End time: %s\n', datestr(now, 'dd/mm/yyyy HH:MM:SS'));
    fprintf(log_file, '- Files timestamp: %s\n', run_timestamp);
    
    fprintf('\nAll files saved with timestamp: %s\n', run_timestamp);
    fprintf('Log file: %s\n', log_filename);
    
    % Close log file
    fclose(log_file);
end

% ========================================================================
% Auxiliary Functions (Same as before but more streamlined)
% ========================================================================

function x = x_coord(a, DeltaX)
    x = a * DeltaX;
end

function y = y_coord(a, DeltaX, GrossStep, Y_terrain)
    pos = a * DeltaX;
    section = floor(pos / GrossStep) + 1;
    
    if section >= length(Y_terrain)
        y = Y_terrain(end);
        return;
    end
    
    dx = pos - (section-1)*GrossStep;
    prop = dx / GrossStep;
    y = Y_terrain(section) + prop*(Y_terrain(section+1)-Y_terrain(section));
end

function R = distance_source_p(p, Xsource, Ysource, X_terrain, Y_terrain, DeltaX, GrossStep)
    xp = x_coord(p, DeltaX);
    yp = y_coord(p, DeltaX, GrossStep, Y_terrain);
    R = sqrt((Xsource - xp)^2 + (Ysource - yp)^2);
end

function R = distance_source_obs(p, Xsource, Ysource, X_terrain, Y_terrain, DeltaX, GrossStep, H)
    xp = x_coord(p, DeltaX);
    yp = y_coord(p, DeltaX, GrossStep, Y_terrain);
    R = sqrt((Xsource - xp)^2 + (Ysource - yp - H)^2);
end

function R = distance_p_q(p, q, X_terrain, Y_terrain, DeltaX, GrossStep)
    xp = x_coord(p, DeltaX);
    yp = y_coord(p, DeltaX, GrossStep, Y_terrain);
    xq = x_coord(q, DeltaX);
    yq = y_coord(q, DeltaX, GrossStep, Y_terrain);
    R = sqrt((xq - xp)^2 + (yq - yp)^2);
end

function R = distance_surf_obs(n, idx, X_terrain, Y_terrain, DeltaX, GrossStep, H)
    xn = x_coord(n, DeltaX);
    yn = y_coord(n, DeltaX, GrossStep, Y_terrain);
    xobs = x_coord(idx, DeltaX);
    yobs = y_coord(idx, DeltaX, GrossStep, Y_terrain) + H;
    R = sqrt((xobs - xn)^2 + (yobs - yn)^2);
end

function L = segment_length(p, X_terrain, Y_terrain, DeltaX, GrossStep)
    if p == 1
        R = distance_p_q(p, p+1, X_terrain, Y_terrain, DeltaX, GrossStep);
    elseif p == floor((10.0 * 70) / DeltaX)
        R = distance_p_q(p-1, p, X_terrain, Y_terrain, DeltaX, GrossStep);
    else
        R1 = distance_p_q(p-1, p, X_terrain, Y_terrain, DeltaX, GrossStep);
        R2 = distance_p_q(p, p+1, X_terrain, Y_terrain, DeltaX, GrossStep);
        R = (R1 + R2) / 2;
    end
    L = R;
end

function write_results_timestamped(J, Et, X_terrain, Y_terrain, DeltaX, GrossStep, NoLinesubs, ...
                                  Xsource, Ysource, H, timestamp)
    % Surface Current Output
    current_filename = sprintf('../results/SurfaceCurrent_MoM_%s.dat', timestamp);
    fileID = fopen(current_filename, 'w');
    fprintf(fileID, '%% MoM Surface Current Results - Generated: %s\n', datestr(now));
    fprintf(fileID, '%% Distance(m)\tMagnitude(A/m)\tReal(A/m)\tImag(A/m)\tPhase(deg)\n');
    for idx = 1:NoLinesubs
        x_pos = x_coord(idx, DeltaX);
        phase_deg = angle(J(idx)) * 180 / pi;
        fprintf(fileID, '%.6f\t%.6e\t%.6e\t%.6e\t%.3f\n', ...
            x_pos, abs(J(idx)), real(J(idx)), imag(J(idx)), phase_deg);
    end
    fclose(fileID);
    
    % Electric Field Output
    field_filename = sprintf('../results/ElectricField_MoM_%s.dat', timestamp);
    fileID = fopen(field_filename, 'w');
    fprintf(fileID, '%% MoM Electric Field Results - Generated: %s\n', datestr(now));
    fprintf(fileID, '%% Distance(m)\tField_dB\tMagnitude(V/m)\tReal(V/m)\tImag(V/m)\tPhase(deg)\n');
    for idx = 1:NoLinesubs
        x_pos = x_coord(idx, DeltaX);
        R_obs = distance_source_obs(idx, Xsource, Ysource, X_terrain, Y_terrain, DeltaX, GrossStep, H);
        dB_field = 20*log10(abs(Et(idx))/sqrt(R_obs));
        phase_deg = angle(Et(idx)) * 180 / pi;
        fprintf(fileID, '%.6f\t%.6f\t%.6e\t%.6e\t%.6e\t%.3f\n', ...
            x_pos, dB_field, abs(Et(idx)), real(Et(idx)), imag(Et(idx)), phase_deg);
    end
    fclose(fileID);
    
    % Quick CSV Summary
    csv_filename = sprintf('../results/Results_Summary_%s.csv', timestamp);
    fileID = fopen(csv_filename, 'w');
    fprintf(fileID, 'Distance_m,Current_Mag,Current_Phase_deg,Field_Mag,Field_dB\n');
    for idx = 1:NoLinesubs
        x_pos = x_coord(idx, DeltaX);
        current_phase = angle(J(idx)) * 180 / pi;
        R_obs = distance_source_obs(idx, Xsource, Ysource, X_terrain, Y_terrain, DeltaX, GrossStep, H);
        dB_field = 20*log10(abs(Et(idx))/sqrt(R_obs));
        fprintf(fileID, '%.6f,%.6e,%.3f,%.6e,%.3f\n', ...
            x_pos, abs(J(idx)), current_phase, abs(Et(idx)), dB_field);
    end
    fclose(fileID);
end

function plot_results_timestamped(J, Et, X_terrain, Y_terrain, DeltaX, GrossStep, NoLinesubs, ...
                                 Xsource, Ysource, H, timestamp)
    % Generate coordinates
    x_coords = arrayfun(@(idx) x_coord(idx, DeltaX), 1:NoLinesubs);
    
    % Main results plot
    fig1 = figure('Units','normalized','Position',[0.1 0.1 0.85 0.75]);
    
    % Current magnitude
    subplot(2,3,1);
    plot(x_coords, abs(J), 'b-', 'LineWidth', 2);
    xlabel('Distance (m)');
    ylabel('|J| (A/m)');
    title('Current Magnitude');
    grid on;
    
    % Current phase
    subplot(2,3,2);
    plot(x_coords, angle(J)*180/pi, 'r-', 'LineWidth', 2);
    xlabel('Distance (m)');
    ylabel('Phase (deg)');
    title('Current Phase');
    grid on;
    
    % Field magnitude
    subplot(2,3,3);
    plot(x_coords, abs(Et), 'g-', 'LineWidth', 2);
    xlabel('Distance (m)');
    ylabel('|E| (V/m)');
    title('Field Magnitude');
    grid on;
    
    % Field in dB
    subplot(2,3,4);
    R_obs = arrayfun(@(idx) distance_source_obs(idx, Xsource, Ysource, X_terrain, Y_terrain, DeltaX, GrossStep, H), 1:NoLinesubs);
    dB_fields = 20*log10(abs(Et)./sqrt(R_obs));
    plot(x_coords, dB_fields, 'c-', 'LineWidth', 2);
    xlabel('Distance (m)');
    ylabel('Field (dB)');
    title('Field in dB');
    grid on;
    
    % Field phase
    subplot(2,3,5);
    plot(x_coords, angle(Et)*180/pi, 'm-', 'LineWidth', 2);
    xlabel('Distance (m)');
    ylabel('Phase (deg)');
    title('Field Phase');
    grid on;
    
    % Terrain setup
    subplot(2,3,6);
    y_terrain_interp = arrayfun(@(idx) y_coord(idx, DeltaX, GrossStep, Y_terrain), 1:NoLinesubs);
    plot(x_coords, y_terrain_interp, 'k-', 'LineWidth', 2);
    hold on;
    plot(Xsource, Ysource, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    plot(x_coords, y_terrain_interp + H, 'b--', 'LineWidth', 1);
    xlabel('Distance (m)');
    ylabel('Height (m)');
    title('Setup');
    legend('Terrain', 'Source', 'Obs.', 'Location', 'best');
    grid on;
    
    sgtitle(sprintf('MoM Solution Results - %s', timestamp), 'FontSize', 14, 'FontWeight', 'bold');
    
    % Save plots
    fig_name = sprintf('../results/MoM_Results_%s', timestamp);
    saveas(fig1, [fig_name '.fig']);
    saveas(fig1, [fig_name '.png']);
    close(fig1);
end