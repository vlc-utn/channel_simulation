%% LOS Channel Single Source
% Simulate LOS channel gain for a single source.
% This code replicates the results from:
%   Visible Light Communication: A Comprehensive Theory and Applications
%   with MATLAB®. Page 52.

clc; clear; close all;

%% Room parameters
lx = 6; ly = 6; lz = 3;             % Dimensions of the Room Environment [m]

%% Tx parameters
half_angle = 70;                    % Semi angle of the LED at half power illumination [degree] (I(half_angle) = 1/2 * I(0), from the Lambertian distribution)
Pt = 20 * 3600e-3;                  % Transmitted optical power from the LED [W]
r_s = [-lx/4, -ly/4, 0];            % Position of LED [m] (z=0 is the roof, z=lz is the floor, x=-lx/2 is the left wall)
n_s = [0, 0, 1];                    % Orientation of the source (z=1 is looking down, z=-1 is looking up)

%% Rx parameters
area = 0.001;                       % Area of the Photodiode [m]
Ts = 1;                             % Gain of the Optical Filter
n = 1.5;                            % Refractive Index of the Lens
FOV = 70;                           % Field of View of the Photodiode
Responsivity = 1;                   % Responsivity [A/W]

% Position of receiver
x = linspace(-lx/2, lx/2, lx*5);    % Points to evaluate in the "X" axis.
y = linspace(-ly/2, ly/2, ly*5);    % Points to evaluate in the "Y" axis.
z = 2.25;                           % Points to evaluate in the "Z" axis.
n_r = [0, 0, -1];                   % Orientation of the receiver (z=1 is looking down, z=-1 is looking up)

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

%% Calculations

% Intermediate variables
m = -log(2)/log(cosd(half_angle));  % Lamberts Mode Number for Tx
g = (n^2)/(sind(FOV).^2);           % Gain of the optical concentrator

% Normalize orientation of senders
for i=1:1:height(n_s)
    n_s(i,:) = n_s(i,:) ./ norm(n_s(i,:));
end

% Normalize orientation of receiver
n_r = n_r ./ norm(n_r);

% 3D coordinates
[XR, YR, ZR] = meshgrid(x, y, z);   % Obtain all possible points in the (X,Y,Z) space
r_r = [XR(:), YR(:), ZR(:)];        % Vectorize. Position of receiver as a 3D vector.

% Pre-allocate vectors
distance = ones(1, length(r_r));
cos_emitter = ones(size(distance));
cos_receiver = ones(size(distance));

% Vector operations
for i=1:1:length(r_r)
    distance(i) = norm(r_s - r_r(i,:));
    cos_emitter(i) = dot(n_s, (r_r(i,:) - r_s) ./ distance(i));
    cos_receiver(i) = dot(n_r, (r_s - r_r(i,:)) ./ distance(i));
end
cos_emitter(cos_emitter < 0) = 0;
cos_receiver(acosd(cos_receiver) > FOV) = 0;
cos_receiver(cos_receiver < 0) = 0;

% Revert to meshgrid coordinates
distance = reshape(distance, size(XR));
cos_emitter = reshape(cos_emitter, size(XR));
cos_receiver = reshape(cos_receiver, size(XR));

% LOS channel response
H_LOS = ( (m+1) / (2*pi) ) .* cos_emitter.^m .* ...
    ( area .* cos_receiver ./ (distance.^2) ) .* Ts .* g;

% Received power
P_received = Pt .* H_LOS * Responsivity;
P_received_dBm = 10*log10(P_received/1e-3);

%% Figure
surfc(x, y, P_received_dBm);
title('Received Power in Indoor - VLC System corresponding to the LOS path');
xlabel('x in m');
ylabel('y in m');
zlabel('Received Power in dBm');
axis([-lx/2, lx/2, -ly/2, ly/2, min(min(P_received_dBm)), max(max(P_received_dBm))]);