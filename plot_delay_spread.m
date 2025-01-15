function plot_delay_spread(x_rx, y_rx, lx, ly, mean_delay, Drms, ...
    Drms_nlos, Bc, Bc_nlos, printPDF, fsize)
%PLOT_DELAY_SPREAD Plot delay spread

    %%% Mean delay total channel
    if(printPDF)
        figa = figure(windowState="maximized");
    else
        figure(NumberTitle="off", Name="Delay Spread");
        subplot(2,2,1);
    end
    surf(x_rx, y_rx, mean_delay/1e-9);
    xlabel('x [m]', Interpreter='latex');
    ylabel('y [m]', Interpreter='latex');
    zlabel("Dispersi\'on de retraso medio $\bar{\tau}$ [ns]", Interpreter='latex');
    axis([-lx/2, lx/2, -ly/2, ly/2, min(min(mean_delay/1e-9)), max(max(mean_delay/1e-9))]);
    colorbar;
    ax = gca;
    ax.FontSize = fsize;
    
    %%% RMS Spread total channel
    if(printPDF)
        figb = figure(windowState="maximized");
    else
        subplot(2,2,2);
    end
    surf(x_rx, y_rx, Drms/1e-9);
    xlabel('x [m]', Interpreter='latex');
    ylabel('y [m]', Interpreter='latex');
    zlabel("Dispersi\'on RMS $\tau_{RMS}$ [ns]", Interpreter='latex');
    axis([-lx/2, lx/2, -ly/2, ly/2, min(min(Drms/1e-9)), max(max(Drms/1e-9))]);
    colorbar;
    ax = gca;
    ax.FontSize = fsize;
    
    %%% Coherence bandwidth total channel
    if(printPDF)
        figd = figure(windowState="maximized");
    else
        subplot(2,2,3);
    end
    surf(x_rx, y_rx, Bc/1e6);
    xlabel('x [m]', Interpreter='latex');
    ylabel('y [m]', Interpreter='latex');
    zlabel('$B_{c}$ [MHz]', Interpreter='latex');
    axis([-lx/2, lx/2, -ly/2, ly/2, min(min(Bc/1e6)), max(max(Bc/1e6))]);
    colorbar;
    ax = gca;
    ax.FontSize = fsize;
    
    %%% Coherence bandwidth NLOS channel
    if(printPDF)
        fige = figure(windowState="maximized");
    else
        subplot(2,2,4);
    end
    surf(x_rx, y_rx, Bc_nlos/1e6);
    xlabel('x [m]', Interpreter='latex');
    ylabel('y [m]', Interpreter='latex');
    zlabel('$B_{c}$ NLOS [MHz]', Interpreter='latex');
    axis([-lx/2, lx/2, -ly/2, ly/2, min(min(Bc_nlos/1e6)), max(max(Bc_nlos/1e6))]);
    colorbar;
    ax = gca;
    ax.FontSize = fsize;

    if(printPDF)
        exportgraphics(figa, 'images/mean_delay.pdf', ContentType='vector')
        exportgraphics(figb, 'images/rms_spread.pdf', ContentType='vector')
        exportgraphics(figd, 'images/coherence_bandwidth.pdf', ContentType='vector')
        exportgraphics(fige, 'images/coherence_bandwidth_nlos.pdf', ContentType='vector')
    end

end

