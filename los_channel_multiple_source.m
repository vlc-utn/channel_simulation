% TODO this file is a copy paste from the book, it should use the same
% syntax from the single source example.

%% LOS Channel Multiple Source
% Simulate LOS channel gain for multiple sources
% This code replicates the results from:
%   Visible Light Communication: A Comprehensive Theory and Applications
%   with MATLAB®. Page 54.

clc; clear; close all ;

%% Parametros del ambiente
lx=6; ly=6; lz=3; % Dimensions of the Room Environment 6 X 6 X 3 m ^{3}
h = 2.25; % Distance between the transmitter and the receiver plane
[XT, YT] = meshgrid ([-lx/4, lx/4], [-ly/4, ly/4]); % Position of LED

x = linspace(-lx/2, lx/2, lx*5);
y = linspace(-ly/2, ly/2, ly*5);
[XR, YR] = meshgrid(x, y ); % Illustrates the receiver

%% Parámetros del transmisor
theta =60; % Semi angle of the LED at half power illumination
m = -log10(2) / log10(cosd(theta)); % Lamberts Mode Number
Power_LED = 20; % Transmitted Power of the LED
No_LED = 3600; % Total nummber of LEDs
Total_Power = No_LED * Power_LED; % Total amount of Power contributed by all the LEDs plane grid

D = sqrt((XR - XT (1 ,1)).^2 + (YR - YT(1 ,1)).^2 + h^2); % Distance
receiver_angle = acosd(h ./ D);

%% Parámetros del receptor
APD = 0.001; % Area of the Photodiode
rho = 0.8; % Reflectance Parameter
Ts = 1; % Gain of the Optical Filter
Refractive_Index = 1.5; % Refractive Index of the Lens
FOV_PD = 70; % Field of View of the Photodiode
Concentrator_Gain = (Refractive_Index^2) / (sind(FOV_PD).^2); % Gain of the optical concentrator

%% Cálculos
H_A1 = (m+1) * APD .* cosd(receiver_angle).^m ./ (2*pi .* D.^2) ; % Channel DC gain corresponding to the LOS path

P_received = Total_Power .* H_A1 .* Ts .* Concentrator_Gain ; % Total amount of Received power corresponding to single LED source 1
P_received(find(abs(receiver_angle) > FOV_PD)) = 0; % Si el ángulo receptor es mayor que el FOV, la potencia es cero.

P_rec_A2 = fliplr ( P_received ); % Computation of Received power from source 2 , due to symmetry separate calculations are not required
P_rec_A3 = flipud ( P_received );
P_rec_A4 = fliplr ( P_rec_A3 );

P_received_total = P_received + P_rec_A2 + P_rec_A3 + P_rec_A4 ; % Total received power from both LED sources
P_received_dBm =10*log10 ( P_received_total );

%% Figure
surfc (x ,y , P_received_dBm ) ; hold on
title ('Received Power in Indoor visible light communication ')
xlabel ( 'x in m ')
ylabel ( 'y in m ')
zlabel ( ' Received Power in dBm ')
% Computation of SNR = Signal Power / Noise Power
noisebandwidth_factor = 0.562; % Noise Bandwidth Factor
dataRate =512000; % Data Rate
Bn = noisebandwidth_factor * dataRate ; % Total Bandwidth
q =1.6*10e-19; % Charge of the electron
R =1; % Responsivity of the Photodiode
Sigma_shot =2* q * R *( P_received_total ) * Bn ; % Computation of Shot Noise
Ba =4.5e6 ;
Amplifier_current =0.01;
Sigma_Amplifier = Amplifier_current ^2* Ba ; % Computation of receiver circuitry noise
Sigma_Total = Sigma_shot + Sigma_Amplifier ; % Total Noise
SNR =( R * P_received_total ) ^2./ Sigma_Total ; % Determination of Signal to Noise Ratio ( SNR )
SNRdb =10* log10 ( SNR ) ; % SNR in dB
figure (2)
surfc (x ,y , SNRdb )
xlabel ( 'x in m ') ;
ylabel ( 'y in m ') ;
zlabel ( ' SNR Distribution in the room in dB ') ;