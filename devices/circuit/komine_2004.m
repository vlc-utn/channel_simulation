%% Receiver circuit from (Komine, 2004)

q = 1.60217663e-19;             % [Coulomb] Electron charge.
Kb = physconst("Boltzmann");    % [J/K] Boltzmann constant.
T = 295;                        % [K] Noise temperature.

I_bg = 5100e-6;                  % [A] Background photocurrent produced by ambient ilumination (Moreira, 1997)"

dia_df = 0;                     % [A^2/Hz] Shunt noise current at the input from the amplifier.
dea_df = 4*Kb*T*1.5/(30e-3);    % [V^2/Hz] Series noise voltage at the input from the amplifier.

Rin = 10;                       % [Ohm] Input impedance
Ca = 0;                         % [Farad] Input capacitance of the circuit (not the photodiode)

B = 100e6;                      % [Hz] Bandwidth
I2 = 0.562;                     % Noise-bandwidth factor.
I3 = 0.0868;                    % Noise-badwidth factor.