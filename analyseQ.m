%{
    Function for analysing Q matrix of a one-step FDS update
    For an FDS (in matrix form) of the form
        
    Au^{n+1} = Bu^n + Cu^{n-1}

    This can be written in one-step form according to

    --       --       --                     --  --       --
    | u^{n+1} |       |   A^{-1}B    A^{-1}C  |  |   u^n   | 
    |         |   =   |                       |  |         |
    |   u^n   |       |      I          0     |  | u^{n-1} | 
    --       --       --                     --  --       -- 

      w^{n+1}     =                Q                 w^n

    
%}
function [f, sigma, z] = analyseQ (Q, k, varargin)

    if length(varargin) == 0
        holdon = false;
        color = 'b';
    end
    if length(varargin) > 0
        holdon = varargin{1};
        if holdon
            color = 'r';
        else
            color = 'b';
        end
    end    
    if length(varargin) > 1
        plotting = varargin{2};
    else 
        plotting = false;
    end
    
    scatterPlot = false;
    
    z = eig(Q);         % Solutions to eigenvalue problem
    s = log(z)/k;       % z = e^{sk} --> s = ln(z)/k

    % As s = j*omega + sigma, we can obtain the (angular) eigenfrequencies
    % by taking the imaginary part of s and the damping coefficients by
    % taking the real part of s. 
    
    % As all coefficients of Q are real, s comes in complex conjugates. We 
    % only want to save the non-negative imaginary values as they
    % correspond to non-negative frequencies.
    
    s = s(imag(s)>0);
    
    [omega, order] = sort(imag(s)); % also save the order of the frequencies to be used for sigma    
    f = omega/(2*pi);      % convert to frequency in Hz
    
    sigma = real(s(order)); % get the damping coefficients 
    
    %% Plot results
    if plotting
        % Eigenfrequencies
        subplot(211)
        if holdon
            hold on;
        else
            hold off;
        end
        if scatterPlot
            scatter(1:length(f), f, color, 'Linewidth', 2)
        else
            plot(1:length(f), f, color, 'Linewidth', 2)
        end

        %figure settings
        title("Eigenfrequencies")
        grid on;
        xlabel("$p$", 'interpreter', 'latex');
        ylabel("$f$ (in Hz)", 'interpreter', 'latex');
        set(gca, 'Fontsize', 16, 'Linewidth', 2)


        % Damping
        subplot(212)
        if holdon
            hold on;
        else
            hold off;
        end

        if scatterPlot
            scatter(1:length(sigma), sigma, color, 'Linewidth', 2)
        else
            plot(1:length(sigma), sigma, color, 'Linewidth', 2)
        end
         %figure settings
        title("Damping")
        grid on;
        xlabel("$p$", 'interpreter', 'latex');
        ylabel("$\sigma$", 'interpreter', 'latex');
        set(gca, 'Fontsize', 16, 'Linewidth', 2)
        set(gcf, 'color', 'w')
    end
end
