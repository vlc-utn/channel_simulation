%% CREE RED LED
% https://ar.mouser.com/ProductDetail/Cree-LED/XPERED-L1-0000-00501?qs=y%252BJdrdj3vZo8w6cewjsMLA%3D%3D
% XPERED-L1-0000-00501

% Taken from datasheet "color" pg. 18
half_angle = 70;                % [degree] Semi angle of the LED at half power illumination (I(half_angle) = 1/2 * I(0), from the Lambertian distribution)

% PPF = 1.48, lambda = 625nm -> PPF = Po / (h*c/lambda) [fotones / seg]
PLANCK = 6.626e-34;
MOL = 6.022e23;
C = 3e8;
LAMBDA_LED = 625e-9;
PPF = 1.48;

Pt = PPF * MOL*1e-6 * PLANCK * C / LAMBDA_LED;       % 0.2835                 % [W] Transmitted optical power from the LED.

% Flux (lm) = 56.8
I0 = 56.8; % [lm] Total luminic power.