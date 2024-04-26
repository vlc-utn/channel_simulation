%% Circuit for the receiver
% Assumed gain of 20k

% Figure 6-12, page 12
dia_df = (3e-12)^2;             % [A^2 / Hz] Shunt noise current at the input from the amplifier.

% eN parameter, page 8
dea_df = (17.8e-9)^2;           % [V^2 / Hz] Series noise voltage at the input from the amplifier.

% page 7. CPD
Rt = 20e3;                      % [Ohm] Transimpedance gain
Rin = 350;                      % [Ohm] Input impedance

% Invented value
Ca = 0.1e-9;                    % [Farad] Input capacitance of the circuit (not the photodiode)

B = 50e6;                       % [Hz] Bandwidth