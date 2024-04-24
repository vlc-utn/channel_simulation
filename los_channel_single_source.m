%% LOS Channel Single Source
% Simulate LOS channel gain for a single source.
% This code replicates the results from:
%   Visible Light Communication: A Comprehensive Theory and Applications
%   with MATLAB®. Page 52.

clc; clear; close all;

%% Room parameters
lx = 6; ly = 6; lz = 3;             % [m] Dimensions of the Room Environment 
rho = 0.7;                          % Reflection coefficient of the wall. Value is from EN 12464-1. Section 4.2.2 "Reflectance of surfaces"

%% Tx parameters
half_angle = 70;                    % [degree] Semi angle of the LED at half power illumination (I(half_angle) = 1/2 * I(0), from the Lambertian distribution)
Pt = 1;                             % [W] Transmitted optical power from the LED.
I0 = 900;                           % [lm] Total luminic power.

% Position of LED [m] (z=0 is the roof, z=lz is the floor, x=-lx/2 is the
% left wall). Accepts multiple sources
r_s = [ -lx/4,  -ly/4,  0;
        -lx/4,  ly/4,   0;
        lx/4,   -ly/4,  0;
        lx/4,   ly/4,   0];

% Orientation of the source (z=1 is looking down, z=-1 is looking up).
% Accepts multiple sources.
n_s = [ 0,  0,  1;
        0,  0,  1;
        0,  0,  1;
        0,  0,  1];

%% Rx parameters
area = 0.001;                       % [m] Area of the Photodiode.
Ts = 1;                             % Gain of the Optical Filter.
n = 1.5;                            % Refractive Index of the Lens.
FOV = 70;                           % [degree] Field of View of the Photodiode.
z_rx = 2.25;                        % [m] Position of the receiver in the "Z" axis.
n_r = [0, 0, -1];                   % Orientation of the receiver (z=1 is looking down, z=-1 is looking up).

%% Simulation parameters
% Number of points to evaluate for the simulation
Nx = round(lx*5);
Ny = round(ly*5);
Nz = round(lz*5);
dA = lz*ly / (Ny*Nz);                   % Differential area used for wall bounces. This value is suggested in (Ghassemlooy 2018), page 90
dt = sqrt(dA)/physconst("LightSpeed");  % Delta time for simulation, as proposed in (Barry, 1993).
t_vector = 0:dt:30e-9;                  % Temporal vector for simulation.


%% Variable check
if (height(r_s) ~= height(n_s))
    error("Different number of position and orientations for Tx. (r_s and n_s don't match)");
end

for i=1:1:height(r_s)
    if ((r_s(i,1) < -lx/2) || (r_s(i,1) > lx/2))
        error("Source '%d' is out of bounds in the X axis. Valid range: [%.2f, %.2f]. Current value: %.2f",i, -lx/2, lx/2, r_s(i,1));
    elseif ((r_s(i,2) < -ly/2) || (r_s(i,2) > ly/2))
        error("Source '%d' is out of bounds in the Y axis. Valid range: [%.2f, %.2f]. Current value: %.2f",i, -ly/2, ly/2, r_s(i,2));
    elseif ((r_s(i,3) < 0) || (r_s(i,3) > lz))
        error("Source '%d' is out of bounds in the Z axis. Valid range: [%.2f, %.2f]. Current value: %.2f",i, 0, lz, r_s(i,3));
    end
end

%% Vector manipulation
% Receiver position
x_rx = linspace(-lx/2, lx/2, Nx);           % Points to evaluate in the "X" axis for the receiver.
y_rx = linspace(-ly/2, ly/2, Ny);           % Points to evaluate in the "Y" axis for the receiver.
[XR, YR, ZR] = meshgrid(x_rx, y_rx, z_rx);  % Obtain all possible points in the (X,Y,Z) space for the receiver
r_r = [XR(:), YR(:), ZR(:)];                % Vectorize. Position of receiver as a 3D vector.

% Coordinates for the four walls.
x_wall = linspace(-lx/2, lx/2, Nx);     % Points to evaluate in the "X" axis for the walls.
y_wall = linspace(-ly/2, ly/2, Ny);     % Points to evaluate in the "Y" axis for the walls.
z_wall = linspace(0, lz, Nz);           % Points to evaluate in the "Z" axis for the walls.

n_lw = [0 1 0];     % Orientation left wall.
n_rw = [0 -1 0];    % Orientation right wall.
n_bw = [1 0 0];     % Orientation back wall.
n_fw = [-1 0 0];    % Orientation front wall.

[X_LW, Y_LW, Z_LW] = meshgrid(x_wall, -ly/2, z_wall);
[X_RW, Y_RW, Z_RW] = meshgrid(x_wall, ly/2, z_wall);
[X_BW, Y_BW, Z_BW] = meshgrid(-lx/2, y_wall, z_wall);
[X_FW, Y_FW, Z_FW] = meshgrid(lx/2, y_wall, z_wall);

r_walls = [X_LW(:), Y_LW(:), Z_LW(:);
           X_RW(:), Y_RW(:), Z_RW(:);
           X_BW(:), Y_BW(:), Z_BW(:);
           X_FW(:), Y_FW(:), Z_FW(:)];

n_walls = [repmat(n_lw, numel(X_LW), 1);
           repmat(n_rw, numel(X_RW), 1);
           repmat(n_bw, numel(X_BW), 1);
           repmat(n_fw, numel(X_FW), 1)];

% Normalize orientation of senders
for i=1:1:height(n_s)
    n_s(i,:) = n_s(i,:) ./ norm(n_s(i,:));
end

% Normalize orientation of receiver
n_r = n_r ./ norm(n_r);
n_r = repmat(n_r, height(r_r), 1);

%% Channel calculations

% Intermediate variables
m = -log(2)/log(cosd(half_angle));  % Lamberts Mode Number for Tx
g = (n^2)/(sind(FOV).^2);           % Gain of the optical concentrator
h_s = ones(height(r_s), 1);         % Initial gain of the channel is unity for the transmitter

% LOS channel response
H_LOS = h_channel(r_s, n_s, m, h_s, r_r, n_r, area, FOV) * Ts * g;

P_optical_received = Pt .* H_LOS;
P_optical_received_dBm = 10*log10(P_optical_received/1e-3);
P_optical_received_dBm = reshape(P_optical_received_dBm, size(XR));

% NLOS channel response
H_TX_TO_WALL = h_channel(r_s, n_s, m, h_s, r_walls, n_walls, dA, 90);
H_WALL_TO_RX = h_channel(r_walls, n_walls, 1, H_TX_TO_WALL, r_r, n_r, area, FOV);
H_NLOS = H_WALL_TO_RX * Ts * g * rho;

P_optical_received_NLOS = Pt .* H_NLOS;
P_optical_received_NLOS_dBm = 10*log10(P_optical_received_NLOS/1e-3);
P_optical_received_NLOS_dBm = reshape(P_optical_received_NLOS_dBm, size(XR));

%% Figure
figure();
surfc(x_rx, y_rx, P_optical_received_dBm);
title('Received Optical Power in Indoor - VLC System corresponding to the LOS path');
xlabel('x in m');
ylabel('y in m');
zlabel('Received Optical Power in dBm');
axis([-lx/2, lx/2, -ly/2, ly/2, min(min(P_optical_received_dBm)), max(max(P_optical_received_dBm))]);

figure();
surfc(x_rx, y_rx, P_optical_received_NLOS_dBm);
title('Received Optical Power in Indoor - VLC System corresponding to the NLOS path');
xlabel('x in m');
ylabel('y in m');
zlabel('Received Optical Power in dBm');
axis([-lx/2, lx/2, -ly/2, ly/2, min(min(P_optical_received_NLOS_dBm)), max(max(P_optical_received_NLOS_dBm))]);

