% physical_constants.m
function constants = physical_constants()
    constants.PI = pi;
    constants.Epsilon_0 = 8.854e-12; % Vacuum permittivity (F/m)
    constants.Mu_0 = 4 * pi * 1e-7;  % Vacuum permeability (H/m)
    constants.c = 1 / sqrt(constants.Mu_0 * constants.Epsilon_0); % Speed of light (m/s)
    constants.f = 970e6;             % Frequency (Hz)
    constants.Lambda = constants.c / constants.f; % Wavelength (m)
    constants.DeltaX = constants.Lambda / 4.0;    % Spatial sampling interval (m)
    constants.Omega = 2.0 * constants.PI * constants.f; % Angular frequency (rad/s)
    constants.Beta_0 = constants.Omega * sqrt(constants.Mu_0 * constants.Epsilon_0); % Wave number (rad/m)
end
