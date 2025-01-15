function plot_temporal(t_vector, H_LOS, H_NLOS1, H_NLOS2, H_NLOS3, ...
    f, DFT, DFT_NLOS, printPDF, fsize, lineW)
%PLOT_TEMPORAL Plot temporal response
    %%% Temporal response LOS + NLOS
    if(printPDF)
        figa = figure(windowState="maximized");
    else
        figure(NumberTitle="off", Name="Temporal Response");
        subplot(2,2,1);
    end

    %%% Calculate power
    power_los   = sum(H_LOS.^2)/length(H_LOS);
    power_nlos1 = sum(H_NLOS1.^2)/length(H_LOS);
    power_nlos2 = sum(H_NLOS2.^2)/length(H_LOS);
    power_nlos3 = sum(H_NLOS3.^2)/length(H_LOS);
    power_nlos_total = power_nlos1 + power_nlos2 + power_nlos3;
    power_total = power_los + power_nlos_total;

    plot(t_vector./1e-9, H_LOS, LineWidth=lineW); hold on;
    plot(t_vector./1e-9, H_NLOS1 + H_NLOS2 + H_NLOS3, LineWidth=lineW); hold on;
    xlabel('Tiempo [ns]', Interpreter='latex');
    ylabel("h(t)  [x$10^{-5}$]", Interpreter='latex');
    legend("$h_{LOS}(t)$", "$h_{NLOS}(t)$", Interpreter='latex', Fontsize=fsize);
    ax = gca;
    ax.FontSize = fsize;
    grid on;

    annotation('textarrow', [0.23 0.20], [0.60 0.55], LineWidth=lineW, ...
        String=sprintf("%.4f%%", power_los/power_total*100), FontSize=fsize);
    annotation('textarrow', [0.28 0.25], [0.18 0.13], LineWidth=lineW, ...
        String=sprintf("%.4f%%", power_nlos_total/power_total*100), FontSize=fsize);

    %%% Temporal response NLOS
    if(printPDF)
        figb = figure(windowState="maximized");
    else
        subplot(2,2,2);
    end
    plot(t_vector./1e-9, H_NLOS1, LineWidth=lineW); hold on;
    plot(t_vector./1e-9, H_NLOS2, LineWidth=lineW); hold on;
    plot(t_vector./1e-9, H_NLOS3, LineWidth=lineW);
    xlabel('Tiempo [ns]', Interpreter='latex');
    ylabel("h(t) [x$10^{-8}$]", Interpreter='latex');
    legend("$h_1(t)$", "$h_2(t)$", "$h_3(t)$", Interpreter='latex', Fontsize=fsize);
    ax = gca;
    ax.FontSize = fsize;
    grid on;

    annotation('textarrow', [0.29 0.26], [0.60 0.55], LineWidth=lineW, ...
        String=sprintf("%.2f%%", power_nlos1/power_nlos_total*100), FontSize=fsize);
    annotation('textarrow', [0.43 0.40], [0.36 0.31], LineWidth=lineW, ...
        String=sprintf("%.2f%%", power_nlos2/power_nlos_total*100), FontSize=fsize);
    annotation('textarrow', [0.55 0.52], [0.22 0.17], LineWidth=lineW, ...
        String=sprintf("%.2f%%", power_nlos3/power_nlos_total*100), FontSize=fsize);
    
    %%% DFT LOS
    if(printPDF)
        figc = figure(windowState="maximized");
    else
        subplot(2,2,3);
    end
    plot(f/1e6, 10*log10(DFT), LineWidth=lineW)
    xlabel("Frecuencia [MHz]", Interpreter='latex');
    ylabel("$|H(f)| [\frac{dB}{Hz}]$", Interpreter='latex');
    grid on;
    ax = gca;
    ax.FontSize = fsize;
    xlim([min(f/1e6), 700]);     % Hardcapped at 700MHz
    
    %%% DFT NLOS
    if(printPDF)
        figd = figure(windowState="maximized");
    else
        subplot(2,2,4);
    end
    plot(f/1e6, 10*log10(DFT_NLOS), LineWidth=lineW);
    xlabel("Frecuencia [MHz]", Interpreter='latex');
    ylabel("$|H_{NLOS}(f)| [\frac{dB}{Hz}]$", Interpreter='latex');
    grid on;
    ax = gca;
    ax.FontSize = fsize;
    xlim([min(f/1e6), 700]); 

    if(printPDF)
        exportgraphics(figa, 'images/temporal_response_LOS.pdf', ContentType='vector');
        exportgraphics(figb, 'images/temporal_response_NLOS.pdf', ContentType='vector');
        exportgraphics(figc, 'images/dft_LOS.pdf', ContentType='vector');
        exportgraphics(figd, 'images/dft_NLOS.pdf', ContentType='vector');
    end
end

