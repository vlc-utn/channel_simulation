% TODO. This file is a copy from the book:
%   Visible Light Communication: A Comprehensive Theory and Applications
%   with MATLAB®. Page 94.

clc; clear; close all;
m =512; % Total number of OFDM symbols
N =1024; % Length of each OFDM symbol
M =4; % Size of the Constellation ( M can be 4 , 8 , 16 , 32 , 64 , 128 , 256)
Ncp =256; % Length of the Cyclic Prefix

% DCO - OFDM Transmitter
Data = randi([0 M-1] , m , N ) ; % Generation of Random bits Matrix of size m by N
DataMod = qammod( Data , M ) ; % Performing Data Modulation
DataMod_serialtoparallel = DataMod.'; % Performing Serial to Parallel Conversion
datamat = DataMod_serialtoparallel; % Assigning the total data to a variable called datamat

% Computation of Hermitian Symmetry Criteria
datamat (1 ,:) =0; % Assigning the First subcarrier to Zero
datamat (513 ,:) =0; % Assigning the Middle Subcarrier to Zero
datamat (514:1024 ,:) = flipud ( conj ( datamat (2:512 ,:) ) ) ; % Illustrating that only half of the subcarriers are exploited for data transmission as the remaining half are flipped complex conjugate versions of the previous ones.
d_ifft = ifft (( datamat ) ) ; % Computation of IFFT operation
d_ifft_paralleltoserial = d_ifft.'; % Parallel to Serial Conversion
CP_part = d_ifft_paralleltoserial(: , end - Ncp +1: end ) ; % Addition of Cyclic Prefix
DCOOFDM_CP =[ CP_part d_ifft_paralleltoserial ]; % Transmissin of DCO - OFDM signal
bdc =7; % DC BIAS
clip = sqrt ((10.^( bdc /10) ) -1) ; % clipping factor k
bdcc = clip * sqrt ( DCOOFDM_CP .* DCOOFDM_CP ) ; % Computation of DC bias
DCOOFDM_BIAS = DCOOFDM_CP + bdcc ; % Addition of DC bias to the cyclic prefix added signal
count =0;
snr_vector =0:1:80; % size of signal to noise ratio (SNR ) vector
for snr = snr_vector
    SNR = snr + 10* log10 ( log2 ( M ) ) ;
    count = count +1 ;
    DCOOFDM_with_chann = awgn(DCOOFDM_BIAS, SNR, 'measured') ; % Addition of AWGN

    % Receiver of DCO - OFDM
    DCOOFDM_with_chann1 = DCOOFDM_with_chann - bdcc ; %Removal of DC bias
    DCOOFDM_removal_CP = DCOOFDM_with_chann1 (: , Ncp +1: N + Ncp ) ; % Removal of Cyclic Prefix
    DCOOFDM_serialtoparallel = DCOOFDM_removal_CP.'; %Serial to Parallel Conversion
    DCOOFDM_parallel_fft = fft ( DCOOFDM_serialtoparallel) ; % Computation of FFT operation
    DCOOFDM_Demodulation = qamdemod (DCOOFDM_parallel_fft.' , M ) ;
    [~ , s_e1(count) ] = symerr ( Data (: ,2:512) , DCOOFDM_Demodulation (: ,2:512) ) ;
end

% %%%%% With 13 dB Bias
m =512; % Total number of OFDM symbols
N =1024; % Length of each OFDM symbol
M =4; % Size of the Constellation ( M can be 4 , 8 , 16 , 32 , 64 , 128 , 256)
Ncp =256 ; % Length of the Cyclic Prefix

% DCO - OFDM Transmitter
Data1 = randi ([0 M-1] , m , N ) ; % Generation of Random bits Matrix of size m by N
DataMod1 = qammod(Data1, M ) ; % Performing Data Modulation
DataMod_serialtoparallel1 = DataMod1.'; % Performing Serial to Parallel Conversion
datamat1 = DataMod_serialtoparallel1; % Assigning the total data to a variable called datamat

% Computation of Hermitian Symmetry Criteria
datamat1 (1 ,:) =0; % Assigning the First subcarrier to Zero
datamat1 (513 ,:) =0; % Assigning the Middle Subcarrier to Zero
datamat1 (514:1024 ,:) = flipud ( conj ( datamat1 (2:512 ,:) ) ) ; % Illustrating that only half of the subcarriers are exploited for data transmission as the remaining half are flipped complex conjugate versions of the previous ones .
d_ifft1 = ifft (( datamat1 ) ) ; % Computation of IFFT operation
d_ifft_paralleltoserial1 = d_ifft1.'; % Parallel to Serial Conversion

CP_part1 = d_ifft_paralleltoserial1 (: , end - Ncp +1: end ) ; % Addition of Cyclic Prefix
DCOOFDM_CP1 =[ CP_part1 d_ifft_paralleltoserial1 ]; % Transmissin of DCO - OFDM signal
bdc1 =13; % DC BIAS
clip1 = sqrt ((10.^( bdc1 /10) ) -1) ; % clipping factor k
bdcc1 = clip1 * sqrt ( DCOOFDM_CP1 .* DCOOFDM_CP1 ) ; % Computation of DC bias
DCOOFDM_BIAS1 = DCOOFDM_CP1 + bdcc1 ; % Addition of DC bias to the cyclic prefix added signal
count =0;
snr_vector =0:1:80; % size of signal to noise ratio (SNR ) vector
for snr = snr_vector
    SNR = snr + 10* log10 ( log2 ( M ) ) ;
    count = count +1 ;
    DCOOFDM_with_channel = awgn (DCOOFDM_BIAS1, SNR, 'measured'); % Addition of AWGN

    % Receiver of DCO - OFDM
    DCOOFDM_with_chann2 = DCOOFDM_with_channel - bdcc1 ; % Removal of DC bias
    DCOOFDM_removal_CP2 = DCOOFDM_with_chann2 (: , Ncp +1: N + Ncp ) ; % Removal of Cyclic Prefix
    DCOOFDM_serialtoparallel2 = DCOOFDM_removal_CP2.'; % Serial to Parallel Conversion
    DCOOFDM_parallel_fft2 = fft (DCOOFDM_serialtoparallel2 ) ; % Computation of FFT operation
    DCOOFDM_Demodulation2 = qamdemod (DCOOFDM_parallel_fft2.' , M ) ;
    [~ , s_e2(count) ]= symerr ( Data1 (: ,2:512) , DCOOFDM_Demodulation2 (: ,2:512) ) ;
end

% Plotting the BER curves
semilogy ( snr_vector , s_e1 , 'rd - ' , 'LineWidth' ,2) ; hold on;
semilogy ( snr_vector , s_e2 , 'gd - ' , 'LineWidth' ,2) ;
legend ('FFT - based DCO - OFDM -7 dB of Bias ' , 'FFT -based DCO - OFDM -13 dB of Bias ') ;
axis ([0 30 10^-4 1]) ;
xlabel ( ' SNR in dB ') ;
ylabel ( ' BER ') ;
grid on ;