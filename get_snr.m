function [snr, S, N] = get_snr(P_optical_total, pd)
    % Get information from circuit and photodetector
    run(pd);

    % Intermediate variables
    %Ct = Ca + Cpd;

    % Signal power
    S = (P_optical_total .* responsivity).^2; % [A^2]

    % Shot Noise
    %sigma_shot = 2*q*responsivity*P_optical_total*B + 2*q*I_bg*I2*B;

    % Thermal noise
    %n_shunt = (4*Kb*T / Rin + dia_df) * B * I2;
    %n_series = dea_df * (B*I2 / Rin^2 + (2*pi*Ct)^2 * B^3 * I3);
    %sigma_thermal  = n_shunt + n_series;

    %n_shunt = (4*Kb*T / Rin + dia_df) * B * I2;
    %n_series = dea_df * (B*I2 / Rin^2 + (2*pi*Ct)^2 * B^3 * I3);
    %sigma_thermal  = n_shunt + n_series;

    sigma_thermal = (4*Kb*T / Rin + dia_df) * B;
    sigma_shot = 2*q*(responsivity*(P_optical_total + pd_noise_power) + I_bg)*B;

    % Total noise
    N = sigma_thermal + sigma_shot;

    % SNR
    snr = S ./ N;

    fprintf("Shot noise: %0.2e [W]\n", mean(mean(sigma_shot)));
    fprintf("Thermal noise: %0.2e [W]\n", mean(mean(sigma_thermal)));
    fprintf("Signal Power: %0.2e [W]\n", mean(mean(S)));
    fprintf("SNR %0.2f [dB]\n", 10*log10(max(max(snr))));
end







