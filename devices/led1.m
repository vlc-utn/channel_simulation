%% LED proposed from "Optical wireless communications: system and channel modelling with MATLAB", page 92, program 3:3

half_angle = 70;    % [degree] Semi angle of the LED at half power illumination (I(half_angle) = 1/2 * I(0), from the Lambertian distribution)
Pt = 1;             % [W] Transmitted optical power from the LED.
I0 = 0.73*3600*2*pi/(1+0.6461);          % [lm] Total luminic power. Value from (Komine, 2004)

