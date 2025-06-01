% physical_constants.m
% This function defines fundamental physical constants.
function constants = physical_constants()
    % constants.PI: Value of Pi (unitless)
    % The ratio of a circle's circumference to its diameter.
    constants.PI = pi;

    % constants.Epsilon_0: Vacuum permittivity (F/m)
    % Represents the capability of a vacuum to permit electric fields.
    constants.Epsilon_0 = 8.854e-12; % F/m

    % constants.Mu_0: Vacuum permeability (H/m)
    % Represents the capability of a vacuum to permit magnetic fields.
    constants.Mu_0 = 4 * pi * 1e-7;  % H/m

    % constants.c: Speed of light in vacuum (m/s)
    % The speed at which all massless particles and associated fields travel in a vacuum.
    % Calculated as 1 / sqrt(Mu_0 * Epsilon_0).
    constants.c = 1 / sqrt(constants.Mu_0 * constants.Epsilon_0); % m/s

    % Note: Frequency-dependent parameters like frequency (f), wavelength (Lambda),
    % spatial sampling interval (DeltaX), angular frequency (Omega), and
    % wave number (Beta_0) have been removed from this file.
    % These parameters are typically application-specific and should be
    % calculated within the specific model or script that requires them,
    % using the model's input frequency. For example, models like Hata
    % take frequency as a direct input.
end
