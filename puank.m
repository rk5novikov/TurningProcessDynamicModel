%      clear all;
%      close all;

    %Initialization of the range of the studied frequencies
    PPP_min = 0.5;
    PPP_max = 5.0;
    PPP_k = 46;
    PPP_step = (PPP_max - PPP_min) / (PPP_k - 1);

    % Poincarization for a system with two degrees of freedom at a depth of cut = 8
    % Initialization of used arrays
    ket = 180;
    oboroty = 70;
    kmax = floor(ket * oboroty);
    XX = zeros(PPP_k, kmax);
    YY = zeros(PPP_k, kmax);
    
    % plot configurations
    i = 1;
    figure;
    set(0, 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman');
    xlabel('p');
    ylabel('A_x', 'Rotation', 0);
    hold on;

    % Plot first graph
    for PPP=PPP_min:PPP_step:PPP_max
        [XX(i,:), YY(i,:)] = main_iter_puankaracii_8(PPP);
        for k = floor(0.5 * kmax):(kmax - 1)
           if ((abs(XX(i, k)) > abs(XX(i, k - 1))) && (abs(XX(i, k)) > abs(XX(i, k + 1))))
                plot(PPP, XX(i, k), 'r.', 'Linewidth', 2);
                grid on
           end;
        end;
        i = i + 1;
    end
    hold off;
    
    % Plot second graph
    i = 1;
    figure;
    set(0, 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman');
    xlabel('p');
    ylabel('A_y','Rotation', 0);
    hold on;
    for PPP=PPP_min:PPP_step:PPP_max
        for k = ( floor( 0.5 * kmax ) ):(kmax - 1)
           if ((abs(YY(i, k)) > abs(YY(i, k - 1))) && (abs(YY(i, k)) > abs(YY(i, k + 1))))
                plot(PPP, YY(i, k), 'r.', 'Linewidth', 2);
                grid on
           end;
        end;
        i = i + 1;
    end
    hold off;

    % Poincarization for a system with two degrees of freedom at a depth of 16
    % Initialization of used arrays
    XX = zeros(PPP_k, kmax);
    YY = zeros(PPP_k, kmax);
    % plot configurations
    i = 1;
    figure;
    set(0, 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman');
    xlabel('p');
    ylabel('A_x', 'Rotation', 0);
    hold on;
    % Poincarization for a system with one degree of freedom
    % Initialization of used arrays
    YY = zeros(PPP_k,kmax);
    % plot configurations
    i = 1;
    figure;
    set(0, 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman');
    xlabel('p');
    ylabel('A_y', 'Rotation', 0);
    hold on;
    % Plot third graph
    for PPP=PPP_min:PPP_step:PPP_max
        YY(i,:) = main_iter_puankaracii_DOF1(PPP);
        for k = ( floor( 0.5 * length(YY(i,:)) ) ):(length( YY(i,:) ) -1)
           if ((abs(YY(i, k)) > abs(YY(i, k - 1))) && (abs(YY(i, k)) > abs(YY(i, k + 1))))
                plot(PPP, YY(i, k), 'r.', 'Linewidth', 2);
                grid on
           end;
        end;
        i = i + 1;
    end
    hold off;


    % Poincarization...
    YY1 = zeros(PPP_k, kmax);
    % Plot options
    i = 1;
    figure;
    set(0, 'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Times New Roman');
    xlabel('p');
    ylabel('A_y', 'Rotation', 0);
    hold on;
    %Plotting
    for PPP=PPP_min:PPP_step:PPP_max
        if( ( (PPP >= 2.6) && (PPP <= 2.8) ) || ( (PPP >= 3.5) && (PPP <= 3.7) ) || ( (PPP >= 4.4) && (PPP <= 4.7) ) )
            YY1(i,1:(floor(ket * (floor(0.58 * oboroty))))) = main_iter_puankaracii_DOF1_8(PPP);
        else
            YY1(i,:) = main_iter_puankaracii_DOF1_8(PPP);
        end
        for k = ( floor( 0.5 * length(YY1(i,:)) ) ):(length( YY1(i,:) ) -1)
           if ((abs(YY1(i, k)) > abs(YY1(i, k - 1))) && (abs(YY1(i, k)) > abs(YY1(i, k + 1))))
             plot(PPP, YY1(i, k), 'r.', 'Linewidth', 2);
             grid on
           end;
        end;
        i = i + 1;
    end
    hold off;