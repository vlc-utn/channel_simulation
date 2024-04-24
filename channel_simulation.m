%% Channel Simulation
% Simulate optical channel.
% This code replicates the results from:
%   Visible Light Communication: A Comprehensive Theory and Applications
%   with MATLAB®. Page 52.
% And:
%   "Optical wireless communications: system and channel modelling with
%   MATLAB". Page 90.
clc; clear; close all;

addpath("devices");
room1;
pd1;
led1;

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
m = -log(2)/log(cosd(half_angle));  % Lamberts Mode Number for Tx
g = (n^2)/(sind(FOV).^2);           % Gain of the optical concentrator
h_s = zeros(height(r_s), length(t_vector)); % Initial gain of the channel is unity for the transmitter
h_s(:,1) = 1;

% LOS channel response
H_LOS = h_channel(r_s, n_s, m, h_s, r_r, n_r, area, FOV, t_vector) * Ts * g;

P_optical_los_dbm = reshape( 10*log10(Pt .* sum(H_LOS, 2) / 1e-3) , size(XR) );

% NLOS first order channel response
H_TX_TO_WALL = h_channel(r_s, n_s, m, h_s, r_walls, n_walls, dA, 90, t_vector);
H_WALL_TO_RX = h_channel(r_walls, n_walls, 1, H_TX_TO_WALL, r_r, n_r, area, FOV, t_vector);
H_NLOS1 = H_WALL_TO_RX * Ts * g * rho;

P_optical_nlos_dbm = reshape( 10*log10(Pt .* sum(H_NLOS1, 2) / 1e-3) , size(XR) );

% NLOS second order channel response
H_WALL_2 = h_channel(r_walls, n_walls, 1, H_TX_TO_WALL, r_walls, n_walls, dA, 90, t_vector);
H_WALL_2_TO_RX = h_channel(r_walls, n_walls, 1, H_WALL_2, r_r, n_r, area, FOV, t_vector);
H_NLOS2 = H_WALL_2_TO_RX * Ts * g *rho^2;

P_optical_nlos2_dbm = reshape( 10*log10(Pt .* sum(H_NLOS2, 2) / 1e-3) , size(XR) );

% NLOS third order channel response
H_WALL_3 = h_channel(r_walls, n_walls, 1, H_WALL_2, r_walls, n_walls, dA, 90, t_vector);
H_WALL_3_TO_RX = h_channel(r_walls, n_walls, 1, H_WALL_3, r_r, n_r, area, FOV, t_vector);
H_NLOS3 = H_WALL_3_TO_RX * Ts * g *rho^3;

P_optical_nlos3_dbm = reshape( 10*log10(Pt .* sum(H_NLOS3, 2) / 1e-3) , size(XR) );

% Total channel response
H_NLOS_TOTAL = H_NLOS1 + H_NLOS2 + H_NLOS3;
H = H_LOS + H_NLOS_TOTAL;
P_optical_total = Pt .* sum(H, 2);
P_optical_total_dbm = 10*log10(P_optical_total/1e-3);
P_optical_total_dbm = reshape(P_optical_total_dbm, size(XR));

% Delay
mean_delay = zeros(1, height(r_r));
Drms = zeros(size(mean_delay));
mean_delay_nlos = zeros(size(mean_delay));
Drms_nlos = zeros(size(mean_delay));
mean_delay_los = zeros(size(mean_delay));
Drms_los = zeros(size(mean_delay));

for j=1:1:height(r_r)
    mean_delay(j) = sum( (H(j,:)).^2 .* t_vector ) / sum(H(j,:).^2);
    Drms(j) = sqrt( sum( (t_vector - mean_delay(j)).^2 .* H(j,:).^2 ) / sum(H(j,:).^2) );

    mean_delay_nlos(j) = sum( (H_NLOS_TOTAL(j,:)).^2 .* t_vector ) / sum(H_NLOS_TOTAL(j,:).^2);
    Drms_nlos(j) = sqrt( sum( (t_vector - mean_delay_nlos(j)).^2 .* H_NLOS_TOTAL(j,:).^2 ) / sum(H_NLOS_TOTAL(j,:).^2));

    mean_delay_los(j) = sum( (H_LOS(j,:)).^2 .* t_vector) / sum(H_LOS(j,:).^2);
    Drms_los(j) = sqrt( sum( (t_vector - mean_delay_los(j)).^2 .* H_LOS(j,:).^2) / sum(H_LOS(j,:).^2));
end

mean_delay = reshape(mean_delay, size(XR));
mean_delay = mean_delay / 1e-9;
Drms = reshape(Drms, size(XR));
Drms = Drms / 1e-9;

mean_delay_nlos = reshape(mean_delay_nlos, size(XR));
mean_delay_nlos = mean_delay_nlos /1e-9;
Drms_nlos = reshape(Drms_nlos, size(XR));
Drms_nlos = Drms_nlos / 1e-9;

mean_delay_los = reshape(mean_delay_los, size(XR));
mean_delay_los = mean_delay_los /1e-9;
Drms_los = reshape(Drms_los, size(XR));
Drms_los = Drms_los / 1e-9;

%% Figure
figure();
subplot(2,2,1);
surfc(x_rx, y_rx, P_optical_los_dbm);
title('Received Optical Power in Indoor - VLC System corresponding to the LOS path');
xlabel('x in m');
ylabel('y in m');
zlabel('Received Optical Power in dBm');
axis([-lx/2, lx/2, -ly/2, ly/2, min(min(P_optical_los_dbm)), max(max(P_optical_los_dbm))]);

subplot(2,2,2);
surfc(x_rx, y_rx, P_optical_nlos_dbm);
title('Received Optical Power in Indoor - VLC System corresponding to the NLOS path');
xlabel('x in m');
ylabel('y in m');
zlabel('Received Optical Power in dBm');
axis([-lx/2, lx/2, -ly/2, ly/2, min(min(P_optical_nlos_dbm)), max(max(P_optical_nlos_dbm))]);

subplot(2,2,3);
surfc(x_rx, y_rx, P_optical_nlos2_dbm);
title('Received Optical Power in Indoor - VLC System corresponding to the NLOS path');
xlabel('x in m');
ylabel('y in m');
zlabel('Received Optical Power in dBm');
axis([-lx/2, lx/2, -ly/2, ly/2, min(min(P_optical_nlos2_dbm)), max(max(P_optical_nlos2_dbm))]);

subplot(2,2,4);
surfc(x_rx, y_rx, P_optical_nlos3_dbm);
title('Received Optical Power in Indoor - VLC System corresponding to the NLOS path');
xlabel('x in m');
ylabel('y in m');
zlabel('Received Optical Power in dBm');
axis([-lx/2, lx/2, -ly/2, ly/2, min(min(P_optical_nlos3_dbm)), max(max(P_optical_nlos3_dbm))]);

figure();
surfc(x_rx, y_rx, mean_delay);
title('Mean delay');
xlabel('x in m');
ylabel('y in m');
zlabel('Mean delay spread [ns]');
axis([-lx/2, lx/2, -ly/2, ly/2, min(min(mean_delay)), max(max(mean_delay))]);

figure();
surfc(x_rx, y_rx, Drms);
title('RMS spread');
xlabel('x in m');
ylabel('y in m');
zlabel('RMS spread [ns]');
axis([-lx/2, lx/2, -ly/2, ly/2, min(min(Drms)), max(max(Drms))]);

figure();
surfc(x_rx, y_rx, Drms_nlos);
title('RMS spread NLOS');
xlabel('x in m');
ylabel('y in m');
zlabel('RMS spread NLOS [ns]');
axis([-lx/2, lx/2, -ly/2, ly/2, min(min(Drms_nlos)), max(max(Drms_nlos))]);

figure();
surfc(x_rx, y_rx, Drms_los);
title('RMS spread LOS');
xlabel('x in m');
ylabel('y in m');
zlabel('RMS spread LOS [ns]');
axis([-lx/2, lx/2, -ly/2, ly/2, min(min(Drms_los)), max(max(Drms_los))]);
