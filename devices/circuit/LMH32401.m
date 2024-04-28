%% Circuit for the receiver (Gain 20K)
q = 1.60217663e-19;             % [Coulomb] Electron charge.
Kb = physconst("Boltzmann");    % [J/K] Boltzmann constant.
T = 300;                        % [K] Noise temperature.

I_bg = 10e-6;                   % [A] Background photocurrent produced by ambient ilumination (Moreira, 1997)"

dia_df = (3e-12)^2;             % [A^2/Hz] Shunt noise current at the input from the amplifier.
dea_df = (17.8e-9)^2;           % [V^2/Hz] Series noise voltage at the input from the amplifier.

Rin = 350;                      % [Ohm] Input impedance
Ca = 0;                         % [Farad] Input capacitance of the circuit (not the photodiode)

B = 50e6;                       % [Hz] Bandwidth
I2 = 0.562;                     % Noise-bandwidth factor.
I3 = 0.0868;                    % Noise-badwidth factor.