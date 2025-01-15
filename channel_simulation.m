%% Channel Simulation
% Simulate optical channel.
% This code replicates the results from:
%   Visible Light Communication: A Comprehensive Theory and Applications
%   with MATLAB®. Page 52.
% And:
%   "Optical wireless communications: system and channel modelling with
%   MATLAB". Page 90.
clc; clear; close all;

%% Input arguments
room = "devices/room/komine_2004";        % Select room environment
led = "devices/leds/komine_2004";                   % Select LED device
pd = "devices/pd/komine_2004";                        % Select photodetector
%room = "devices/room/single_led.m";
%led = "devices/leds/cree_red";
%pd = "devices/pd/sfh_203";

pointToPlot = [0 0 2.15];                              % Point in (x,y,z) coordinates to plot temporal and frequency response

printPDF = false;                                   % Prints PDF files of the graphs
fsize = 32;     % Font Size for axes labels
lineW = 2;      % LineWidth for 2D Graphs

%% Variable check
run(room);
run(pd);
run(led);

% Reset to default fontsize and lineWidth
if(~printPDF)
    fsize = 12;
    lineW = 1;
end

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

x_rx = linspace(-lx/2, lx/2, Nx);           % Points to evaluate in the "X" axis for the receiver.
y_rx = linspace(-ly/2, ly/2, Ny);           % Points to evaluate in the "Y" axis for the receiver.
[XR, YR, ZR] = meshgrid(x_rx, y_rx, z_rx);  % Obtain all possible points in the (X,Y,Z) space for the receiver
r_r = [XR(:), YR(:), ZR(:)];                % Vectorize. Position of receiver as a 3D vector.

% Normalize orientation of senders
for i=1:1:height(n_s)
    n_s(i,:) = n_s(i,:) ./ norm(n_s(i,:));
end

% Normalize orientation of receiver
n_r = n_r ./ norm(n_r);
n_r = repmat(n_r, height(r_r), 1);

%% Channel calculations

% Intermediate variables
m = -log10(2)/log10(cosd(half_angle));  % Lamberts Mode Number for Tx
h_s = zeros(height(r_s), length(t_vector)); % Initial gain of the channel is unity for the transmitter
h_s(:,1) = 1;

% LOS channel response
H_LOS = h_channel(r_s, n_s, m, h_s, r_r, n_r, area, FOV, t_vector) * Ts * g;

P_optical_los = reshape(Pt .* sum(H_LOS, 2), size(XR));
P_optical_los_dbm = 10*log10(P_optical_los / 1e-3);

% NLOS first order channel response
H_TX_TO_WALL = h_channel(r_s, n_s, m, h_s, r_walls, n_walls, dA, 90, t_vector);
H_WALL_TO_RX = h_channel(r_walls, n_walls, 1, H_TX_TO_WALL, r_r, n_r, area, FOV, t_vector);
H_NLOS1 = H_WALL_TO_RX * Ts * g * rho;

P_optical_nlos1 = reshape(Pt .* sum(H_NLOS1, 2), size(XR));
P_optical_nlos1_dbm = 10*log10(P_optical_nlos1 / 1e-3);

% NLOS second order channel response
H_WALL_2 = h_channel(r_walls, n_walls, 1, H_TX_TO_WALL, r_walls, n_walls, dA, 90, t_vector);
H_WALL_2_TO_RX = h_channel(r_walls, n_walls, 1, H_WALL_2, r_r, n_r, area, FOV, t_vector);
H_NLOS2 = H_WALL_2_TO_RX * Ts * g *rho^2;

P_optical_nlos2 = reshape(Pt .* sum(H_NLOS2, 2), size(XR));
P_optical_nlos2_dbm = 10*log10(P_optical_nlos2 / 1e-3);

% NLOS third order channel response
H_WALL_3 = h_channel(r_walls, n_walls, 1, H_WALL_2, r_walls, n_walls, dA, 90, t_vector);
H_WALL_3_TO_RX = h_channel(r_walls, n_walls, 1, H_WALL_3, r_r, n_r, area, FOV, t_vector);
H_NLOS3 = H_WALL_3_TO_RX * Ts * g *rho^3;

P_optical_nlos3 = reshape(Pt .* sum(H_NLOS3, 2), size(XR));
P_optical_nlos3_dbm = 10*log10(P_optical_nlos3 / 1e-3);

% Total channel response
H_NLOS_TOTAL = H_NLOS1 + H_NLOS2 + H_NLOS3;
H = H_LOS + H_NLOS_TOTAL;

P_optical_total = reshape(Pt .* sum(H, 2), size(XR));
P_optical_total_dbm = 10*log10(P_optical_total / 1e-3);

P_optical_nlos_total = reshape(Pt .* sum(H_NLOS_TOTAL, 2), size(XR));
P_optical_nlos_total_dbm = 10*log10(P_optical_nlos_total / 1e-3);

% Iluminance
Iluminance = reshape(get_iluminance(r_s, n_s, m, r_r, n_r, I0), size(XR));

%% Delay calculations
mean_delay = sum( H.^2 .* t_vector, 2) ./ sum(H.^2, 2);
Drms = sqrt(sum((t_vector-mean_delay).^2 .* H.^2, 2) ./ sum(H.^2, 2));
mean_delay = reshape(mean_delay, size(XR));
Drms = reshape(Drms, size(XR));

mean_delay_nlos = sum( H_NLOS_TOTAL.^2 .* t_vector, 2) ./ sum(H_NLOS_TOTAL.^2, 2);
Drms_nlos = sqrt(sum((t_vector-mean_delay_nlos).^2 .* H_NLOS_TOTAL.^2, 2) ./ sum(H_NLOS_TOTAL.^2, 2));
mean_delay_nlos = reshape(mean_delay_nlos, size(XR));
Drms_nlos = reshape(Drms_nlos, size(XR));

Bc = 1./(10*Drms);                % Ancho de banda del canal, según "Optical wireless communications: system and channel modelling with MATLAB, pág 90"
                              % y según "Visible Light Communications, pág 84"
Bc_nlos = 1./(10*Drms_nlos);

% Truncate big edge values
Bc(Bc > 200e6) = 200e6; 

%% SNR calculations
snr = get_snr(P_optical_total, pd);
snr_db = 10*log10(snr);

%% Frequency response
% Select point to be plotted
[~, index] = min(vecnorm(pointToPlot - r_r, 2, 2));

% Calculate PSD
[DFT, f] = pwelch(H(index,:), rectwin(length(H(index,:))), [], 4096, 1/dt, "onesided", "psd");
[DFT_NLOS, ~] = pwelch(H_NLOS_TOTAL(index,:), rectwin(length(H(index,:))), [], 4096, 1/dt, "onesided", "psd");

%% Plotting
plot_optical_power(x_rx, y_rx, lx, ly, Iluminance, P_optical_los_dbm, ...
   P_optical_total_dbm, P_optical_nlos_total_dbm, printPDF, fsize, lineW);

plot_temporal(t_vector, H_LOS(index,:), H_NLOS1(index,:), H_NLOS2(index,:), H_NLOS3(index,:), ...
        f, DFT, DFT_NLOS, printPDF, fsize, lineW);

plot_delay_spread(x_rx, y_rx, lx, ly, mean_delay, Drms, ...
   Drms_nlos, Bc, Bc_nlos, printPDF, fsize);

%% Figure SNR
figure(NumberTitle="off", Name="SnR");
surfc(x_rx, y_rx, snr_db);
title('SnR');
xlabel('x [m]');
ylabel('y [m]');
zlabel('SnR [dB]');
axis([-lx/2, lx/2, -ly/2, ly/2, min(min(snr_db)), max(max(snr_db))]);
