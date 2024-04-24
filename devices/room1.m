%% Room environment from "Optical wireless communications: system and channel modelling with MATLAB", page 92, program 3:3

lx = 5; ly = 5; lz = 3;     % [m] Dimensions of the Room Environment 
rho = 0.8;                  % Reflection coefficient of the wall. Value is from EN 12464-1. Section 4.2.2 "Reflectance of surfaces"

% Position of LED [m] (z=0 is the roof, z=lz is the floor, x=-lx/2 is the
% left wall). Accepts multiple sources
r_s = [ 0,  0,  0];

% Orientation of the source (z=1 is looking down, z=-1 is looking up).
% Accepts multiple sources.
n_s = [ 0,  0,  1];

z_rx = 3;               % [m] Position of the receiver in the "Z" axis.
n_r = [0, 0, -1];       % Orientation of the receiver (z=1 is looking down, z=-1 is looking up).

% Number of points to evaluate for the simulation
Nx = round(lx*5);
Ny = round(ly*5);
Nz = round(lz*5);

dA = lz*ly / (Ny*Nz);                   % Differential area used for wall bounces. This value is suggested in (Ghassemlooy 2018), page 90
dt = sqrt(dA)/physconst("LightSpeed");  % Delta time for simulation, as proposed in (Barry, 1993).
t_vector = 0:dt:40e-9;                  % Temporal vector for simulation.
