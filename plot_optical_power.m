function plot_optical_power(x_rx, y_rx, lx, ly, Iluminance, P_optical_los_dbm, ...
    P_optical_total_dbm, P_optical_nlos_total_dbm, printPDF, fsize, lineW)
%PLOT_OPTICAL_POWER Plot optical power
    
    %%% Iluminance
    if(printPDF)
        figa = figure(windowState="maximized");
    else
        figure(NumberTitle="off", Name="Optical Power");
        subplot(2,2,1);
    end
    surf(x_rx, y_rx, Iluminance);
    xlabel('x [m]', Interpreter='latex');
    ylabel('y [m]', Interpreter='latex');
    zlabel('Iluminancia [lx]', Interpreter='latex');
    axis([-lx/2, lx/2, -ly/2, ly/2, min(min(Iluminance)), max(max(Iluminance))]);
    colorbar;
    ax = gca;
    ax.FontSize = fsize;

    %%% LOS Optical Power
    if(printPDF)
        figb = figure(windowState="maximized");
    else
        subplot(2,2,2);
    end
    surf(x_rx, y_rx, P_optical_los_dbm);
    xlabel('x [m]', Interpreter='latex');
    ylabel('y [m]', Interpreter='latex');
    zlabel('Potencia luminica LOS [dBm]', Interpreter='latex');
    axis([-lx/2, lx/2, -ly/2, ly/2, min(min(P_optical_los_dbm)), max(max(P_optical_los_dbm))]);
    colorbar;
    ax = gca;
    ax.FontSize = fsize;

    %%% NLOS optical power
    if(printPDF)
        figc = figure(windowState="maximized");
    else
        subplot(2,2,3);
    end
    surf(x_rx, y_rx, P_optical_nlos_total_dbm);
    xlabel('x [m]', Interpreter='latex');
    ylabel('y [m]', Interpreter='latex');
    zlabel('Potencia luminica NLOS [dBm]', Interpreter='latex');
    axis([-lx/2, lx/2, -ly/2, ly/2, min(min(P_optical_nlos_total_dbm)), max(max(P_optical_nlos_total_dbm))]);
    colorbar;
    ax = gca;
    ax.FontSize = fsize;

    %%% Total Optical power
    if(printPDF)
        figf = figure(windowState="maximized");
    else
        subplot(2,2,4);
    end
    surf(x_rx, y_rx, P_optical_total_dbm);
    xlabel('x [m]', Interpreter='latex');
    ylabel('y [m]', Interpreter='latex');
    zlabel('Potencia luminica total [dBm]', Interpreter='latex');
    axis([-lx/2, lx/2, -ly/2, ly/2, min(min(P_optical_total_dbm)), max(max(P_optical_total_dbm))]);
    colorbar;
    ax = gca;
    ax.FontSize = fsize;

    if(printPDF)
        exportgraphics(figa, 'images/illuminance.pdf', ContentType='vector')
        exportgraphics(figb, 'images/los_power.pdf', ContentType='vector')
        exportgraphics(figc, 'images/nlos_power.pdf', ContentType='vector')
        exportgraphics(figf, 'images/total_power.pdf', ContentType='vector')
    end
end

