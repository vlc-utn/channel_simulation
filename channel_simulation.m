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
room = "devices/room/channel_modeling_2018";        % Select room environment
led = "devices/leds/generic_led";                   % Select LED device
pd = "devices/pd/VEMD5510C";                        % Select photodetector
circuit = "devices/circuit/LMH32401";               % Select receiving circuit

%% Variable check
run(room);
run(pd);
run(led);
run(circuit);

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

P_optical_los = reshape(Pt .* sum(H_LOS, 2), size(XR));
P_optical_los_dbm = 10*log10(P_optical_los / 1e-3);

% NLOS first order channel response
H_TX_TO_WALL = h_channel(r_s, n_s, m, h_s, r_walls, n_walls, dA, 90, t_vector);
H_WALL_TO_RX = h_channel(r_walls, n_walls, 1, H_TX_TO_WALL, r_r, n_r, area, FOV, t_vector);
H_NLOS1 = H_WALL_TO_RX * Ts * g * rho;

P_optical_nlos = reshape(Pt .* sum(H_NLOS1, 2), size(XR));
P_optical_nlos_dbm = 10*log10(P_optical_nlos / 1e-3);

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

% Iluminance
Iluminance = reshape(get_iluminance(r_s, n_s, m, r_r, n_r, I0), size(XR));

%% Delay calculations
mean_delay = sum( H.^2 .* t_vector, 2) ./ sum(H.^2, 2);
Drms = sqrt(sum((t_vector-mean_delay).^2 .* H.^2, 2) ./ sum(H.^2, 2));
mean_delay = reshape(mean_delay, size(XR)) / 1e-9;
Drms = reshape(Drms, size(XR)) / 1e-9;

mean_delay_nlos1 = sum( H_NLOS1.^2 .* t_vector, 2) ./ sum(H_NLOS1.^2, 2);
Drms_nlos1 = sqrt(sum((t_vector-mean_delay_nlos1).^2 .* H_NLOS1.^2, 2) ./ sum(H_NLOS1.^2, 2));
mean_delay_nlos1 = reshape(mean_delay_nlos1, size(XR)) / 1e-9;
Drms_nlos1 = reshape(Drms_nlos1, size(XR)) / 1e-9;

mean_delay_nlos2 = sum( (H_NLOS1 + H_NLOS2).^2 .* t_vector, 2) ./ sum((H_NLOS1 + H_NLOS2).^2, 2);
Drms_nlos2 = sqrt(sum((t_vector-mean_delay_nlos2).^2 .* (H_NLOS1 + H_NLOS2).^2, 2) ./ sum((H_NLOS1 + H_NLOS2).^2, 2));
mean_delay_nlos2 = reshape(mean_delay_nlos2, size(XR)) / 1e-9;
Drms_nlos2 = reshape(Drms_nlos2, size(XR)) / 1e-9;

mean_delay_nlos3 = sum( H_NLOS_TOTAL.^2 .* t_vector, 2) ./ sum(H_NLOS_TOTAL.^2, 2);
Drms_nlos3 = sqrt(sum((t_vector-mean_delay_nlos3).^2 .* H_NLOS_TOTAL.^2, 2) ./ sum(H_NLOS_TOTAL.^2, 2));
mean_delay_nlos3 = reshape(mean_delay_nlos3, size(XR)) / 1e-9;
Drms_nlos3 = reshape(Drms_nlos3, size(XR)) / 1e-9;

Bc = 1./(10*Drms*1e-9)./1e6;  % Ancho de banda del canal, según "Optical wireless communications: system and channel modelling with MATLAB, pág 90"
                              % y según "Visible Light Communications, pág 84"

% Trim corner values that tend to infinity
Bc(Bc > 500) = 500;

%% SNR calculations
snr = get_snr(P_optical_total, circuit, pd);
snr_db = 10*log(snr);

%% Figure optical power
figure(NumberTitle="off", Name="Optical Power");
subplot(3,2,1);
surfc(x_rx, y_rx, Iluminance);
title('Iluminance');
xlabel('x [m]');
ylabel('y [m]');
zlabel('I [lx]');
axis([-lx/2, lx/2, -ly/2, ly/2, min(min(Iluminance)), max(max(Iluminance))]);

subplot(3,2,2);
surfc(x_rx, y_rx, P_optical_los_dbm);
title('LOS Optical Power');
xlabel('x [m]');
ylabel('y [m]');
zlabel('Optical Power [dBm]');
axis([-lx/2, lx/2, -ly/2, ly/2, min(min(P_optical_los_dbm)), max(max(P_optical_los_dbm))]);

subplot(3,2,3);
surfc(x_rx, y_rx, P_optical_nlos_dbm);
title('NLOS 1st Reflection');
xlabel('x [m]');
ylabel('y [m]');
zlabel('Optical Power [dBm]');
axis([-lx/2, lx/2, -ly/2, ly/2, min(min(P_optical_nlos_dbm)), max(max(P_optical_nlos_dbm))]);

subplot(3,2,4);
surfc(x_rx, y_rx, P_optical_nlos2_dbm);
title('NLOS 2nd Reflection');
xlabel('x [m]');
ylabel('y [m]');
zlabel('Optical Power [dBm]');
axis([-lx/2, lx/2, -ly/2, ly/2, min(min(P_optical_nlos2_dbm)), max(max(P_optical_nlos2_dbm))]);

subplot(3,2,5);
surfc(x_rx, y_rx, P_optical_nlos3_dbm);
title('NLOS 3rd Reflection');
xlabel('x [m]');
ylabel('y [m]');
zlabel('Optical Power [dBm]');
axis([-lx/2, lx/2, -ly/2, ly/2, min(min(P_optical_nlos3_dbm)), max(max(P_optical_nlos3_dbm))]);

subplot(3,2,6);
surfc(x_rx, y_rx, P_optical_total_dbm);
title('Total Optical Power');
xlabel('x [m]');
ylabel('y [m]');
zlabel('Optical Power [dBm]');
axis([-lx/2, lx/2, -ly/2, ly/2, min(min(P_optical_total_dbm)), max(max(P_optical_total_dbm))]);

%% Figure delay spread
figure(NumberTitle="off", Name="Delay Spread");
subplot(3,2,1);
surfc(x_rx, y_rx, mean_delay);
title('Mean delay');
xlabel('x [m]');
ylabel('y [m]');
zlabel('Mean delay spread [ns]');
axis([-lx/2, lx/2, -ly/2, ly/2, min(min(mean_delay)), max(max(mean_delay))]);

subplot(3,2,2);
surfc(x_rx, y_rx, Drms);
title('Total RMS Spread');
xlabel('x [m]');
ylabel('y [m]');
zlabel('RMS Spread [ns]');
axis([-lx/2, lx/2, -ly/2, ly/2, min(min(Drms)), max(max(Drms))]);

subplot(3,2,3);
surfc(x_rx, y_rx, Drms_nlos1);
title('RMS Spread with NLOS1');
xlabel('x [m]');
ylabel('y [m]');
zlabel('RMS Spread [ns]');
axis([-lx/2, lx/2, -ly/2, ly/2, min(min(Drms_nlos1)), max(max(Drms_nlos1))]);

subplot(3,2,4);
surfc(x_rx, y_rx, Drms_nlos2);
title('RMS Spread with NLOS2');
xlabel('x [m]');
ylabel('y [m]');
zlabel('RMS Spread [ns]');
axis([-lx/2, lx/2, -ly/2, ly/2, min(min(Drms_nlos2)), max(max(Drms_nlos2))]);

subplot(3,2,5);
surfc(x_rx, y_rx, Drms_nlos3);
title('RMS Spread with NLOS3');
xlabel('x [m]');
ylabel('y [m]');
zlabel('RMS Spread [ns]');
axis([-lx/2, lx/2, -ly/2, ly/2, min(min(Drms_nlos3)), max(max(Drms_nlos3))]);

subplot(3,2,6);
surfc(x_rx, y_rx, Bc);
title('Channel bandwidth for flat response');
xlabel('x [m]');
ylabel('y [m]');
zlabel('Bc[MHz]');
axis([-lx/2, lx/2, -ly/2, ly/2, min(min(Bc)), max(max(Bc))]);

%% Figure temporal response
figure(NumberTitle="off", Name="Temporal Response");
point = [0 0 3];                                % Point to be plotted
dist = vecnorm(point - r_r, 2, 2);
[~, index] = min(dist);                         % Selecting point to plot

N0 = 1024*2^(nextpow2(length(H(index,:))));     % Samples for the fft
df = 1/(N0*dt) / 1e6;
freq= df*(-N0/2:N0/2-1);                        % Frequency vector

DFT = abs(fftshift(fft(H(index,:), N0)));
DFT_dB = 10*log10(DFT);

DFT_NLOS = abs(fftshift(fft(H_NLOS_TOTAL(index,:), N0)));
DFT_NLOS_dB = 10*log10(DFT_NLOS);

subplot(2,2,1);
plot(t_vector./1e-9, H_LOS(index,:), t_vector./1e-9, H_NLOS1(index,:), t_vector./1e-9, H_NLOS2(index,:), t_vector./1e-9, H_NLOS3(index,:));
title(sprintf("H-LOS(t) at (%0.2f; %0.2f; %0.2f)", r_r(index,1), r_r(index,2), r_r(index,3)));
xlabel('time [ns]');
ylabel('h(t)');
legend("h-los(t)", "h_nlos1(t)", "h_nlos2(t)", "h_nlos3(t)");
grid on;

subplot(2,2,2);
plot(t_vector./1e-9, H_NLOS1(index,:), t_vector./1e-9, H_NLOS2(index,:), t_vector./1e-9, H_NLOS3(index,:));
title(sprintf("H-NLOS(t) at (%0.2f; %0.2f; %0.2f)", r_r(index,1), r_r(index,2), r_r(index,3)));
xlabel('time [ns]');
ylabel('h(t)');
legend("h-nlos1(t)", "h-nlos2(t)", "h-nlos3(t)");
grid on;

subplot(2,2,3);
plot(freq, DFT_dB)
xlabel("Freq [MHz]");
ylabel("|H(f)| [dB]");
title("DFT of |H(f)|");
grid on;

subplot(2,2,4);
plot(freq, DFT_NLOS_dB);
xlabel("Freq [MHz]");
ylabel("|H-NLOS(f)| [dB]");
title("DFT of |H-NLOS(f)|");
grid on;

%% Figure SNR
figure(NumberTitle="off", Name="SnR");
surfc(x_rx, y_rx, snr_db);
title('SnR');
xlabel('x [m]');
ylabel('y [m]');
zlabel('SnR [dB]');
axis([-lx/2, lx/2, -ly/2, ly/2, min(min(snr_db)), max(max(snr_db))]);

