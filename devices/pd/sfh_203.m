%% PD SFH_203
% https://ar.mouser.com/ProductDetail/ams-OSRAM/SFH-203?qs=nTDll3UaDK5tSfK%252By9U02Q%3D%3D

% 1mm^2
area = 1e-6;            % [m] Area of the Photodiode.

n = 1.5;                % Refractive Index of the Lens.

% Usando el LED CREE con lambda = 625, y a partir del gráfico de pg.6
Ts = 0.65;              % Gain of the Optical Filter

% Half angle is 20°. Very not lambertian
FOV = 20;               % [degree] Field of View of the Photodiode.
responsivity = 0.62;    % [A/W] Responsivity of the photodiode.

Cpd = 11e-12;          % [F] Photodiode capacitance.

I_bg = 40e-6;           % [A] Background photocurrent produced by ambient ilumination (Moreira, 1997)"

g = 1;  % Gain of the optical concentrator

%% Circuit for the receiver
% Obtenido del integrado LMH34400
% https://ar.mouser.com/ProductDetail/Texas-Instruments/LMH34400IDRLR?qs=By6Nw2ByBD2lP2c%252BklzHtA%3D%3D
q = 1.60217663e-19;             % [Coulomb] Electron charge.
Kb = physconst("Boltzmann");    % [J/K] Boltzmann constant.
T = 300;                        % [K] Noise temperature.
B = 50e6;
pd_noise_power = 0.029e-12;
dia_df = (3e-12)^2;
Rin = 100;

%% Circuit for the receiver
% Obtenido del integrado LTC6268-10/LTC6269-10
% https://www.mouser.es/ProductDetail/Analog-Devices/LTC6268IS6-10TRMPBF?qs=oahfZPh6IAK%2Fgwt4Y8oAdg%3D%3D
% q = 1.60217663e-19;             % [Coulomb] Electron charge.
% Kb = physconst("Boltzmann");    % [J/K] Boltzmann constant.
% T = 300;                        % [K] Noise temperature.
% 
% dia_df = (7e-15)^2;             % [A^2/Hz] Shunt noise current at the input from the amplifier.
% dea_df = (4e-9)^2;              % [V^2/Hz] Series noise voltage at the input from the amplifier.
% sigma_pd = (0.029e-12)^2;       % [W^2/Hz] From the datasheet
% 
% Rin = 10e3;                      % [Ohm] Input impedance
% Ca = 0.45e-12;                   % [Farad] Input capacitance of the circuit (not the photodiode)
% 
% B = 50e6;                       % [Hz] Bandwidth
% I2 = 0.562;                     % Noise-bandwidth factor.
% I3 = 0.0868;                    % Noise-badwidth factor.