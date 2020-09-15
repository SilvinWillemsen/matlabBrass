%{
    First order brass divided into two parts for dynamic changes
%}

clear all;
close all;

% drawing variables
drawThings = true;
drawStart = 1;
drawSpeed = 10;
centered = true;

fs = 44100;             % Sample rate (Hz)
k = 1/fs;               % Time step (s)
lengthSound = fs * 2;   % Duration (s)

%% viscothermal effects
T = 26.85;
[c, rho, eta, nu, gamma] = calcThermoDynConstants (T);
cInit = c;

%% Tube variables
h = c * k;              % Grid spacing (m)

Ninit = 60.001;
L = Ninit * h;          % Length
N = floor(L/h);         % Number of points (-)
alf = Ninit - N;
N = floor(L/h);         % Number of points (-)

lambda = c * k / h      % courant number

%% Set cross-sectional geometry
[S, SHalf, SBar] = setTube (N+1);

% Quick note: N is the number of spaces between the points so the number of points is N+1

%% Initialise states
upNext = zeros(ceil(N/2) + 1, 1);
up = zeros(ceil(N/2) + 1, 1);
uvNext = zeros(ceil(N/2), 1); 
uv = zeros(ceil(N/2), 1);

up(floor(length(up) / 4 - 5):floor(length(up)/4) + 5) = hann(11);

wpNext = zeros(floor(N/2) + 1, 1);
wp = zeros(floor(N/2) + 1, 1);
wvNext = zeros(floor(N/2), 1); % the total v vector is one smaller than the total p vector
wv = zeros(floor(N/2), 1);

% Initialise output
out = zeros (lengthSound, 1);
outputPos = floor(4/5 * N);

% Set ranges
upRange = 2:length(up)-1;         % range without boundaries
wpRange = 2:length(wp)-1;

%% Initialise energies
kinEnergyU = zeros (lengthSound, 1);
potEnergyU = zeros (lengthSound, 1);
totEnergyU = zeros (lengthSound, 1);

kinEnergyW = zeros (lengthSound, 1);
potEnergyW = zeros (lengthSound, 1);
totEnergyW = zeros (lengthSound, 1);

totEnergy = zeros (lengthSound, 1);
scaledTotEnergy = zeros (lengthSound, 1);

potScalingU = ones(length(up),1);
potScalingU(1) = 1 / 2;
potScalingU(end) = 1 / 2;

potScalingW = ones(length(wp),1);
potScalingW(1) = 1 / 2;
potScalingW(end) = 1 / 2;

ip = [alf * (alf - 1) * (alf - 2) / -6, ...
                (alf - 1) * (alf + 1) * (alf - 2) / 2, ...
                alf * (alf + 1) * (alf - 2) / -2, ...
                alf * (alf + 1) * (alf - 1) / 6];
         
SBarI = SBar(length(uv)-1:length(uv)+2) .* ip';
uvMph = 0;
wvmh = 0;

changeC = false;
for n = 1:lengthSound
    
    % change wave speed
    if changeC
        c = cInit * (1+ n/fs);%* sin(0.5 * pi * n/fs));
    else
        c = c;
    end
    cSave(n) = c;
    
    % save previous state for comparison later
    NPrev = N;

    % recalculate gridspacing, points lambda^2 and alpha from new wave speed
    h = c*k;
    Ninit = L/h;
    N = floor(L/h);
    Nsave(n) = N;
    hSave(n) = h;
    
    lambda = c * k / h;

    alf = (Ninit - N);
  
    % calculate interpolator
%     if interpol == "cubic"
        ip = [alf * (alf - 1) * (alf - 2) / -6, ...
                (alf - 1) * (alf + 1) * (alf - 2) / 2, ...
                alf * (alf + 1) * (alf - 2) / -2, ...
                alf * (alf + 1) * (alf - 1) / 6];
%     else
%         ip = [0, (1-alf), alf, 0];
%     end

    if abs(N - NPrev) > 1
        disp('too fast')
    end
    
    % add point if N^n > N^{n-1}
    if N > NPrev
        customIp = [-alf * (alf + 1) / ((alf + 2) * (alf + 3)), ...
                    2 * alf * (alf + 1) / ((alf + 1) * (alf + 2)), ...
                    2 * (alf + 1) / ((alf + 2) * (alf + 1)), ...
                    -2 * alf / ((alf + 3) * (alf + 2))];
        if mod(N,2) == 1
%             uvNext = [uvNext; (ip(4) * uvNext(end-1) + ip(3) * uvNext(end) + ip(2) * wvNext(1) + ip(1) * wvNext(2))];
%             uv = [uv; (ip(4) * uv(end-1) + ip(3) * uv(end) + ip(2) * wv(1) + ip(1) * wv(2))];
            uvNext = [uv; uvNextMph];
            uv = [uv; uvMph];
            uvMph =  customIp * [uv(end-1:end); wv(1:2)];% should be properly interpolated
            upNext = [upNext; customIp * [upNext(end-1:end); wpNext(1:2)]];
            up = [up; customIp * [up(end-1:end); wp(1:2)]];
            
        else 
%             wvNext = [(ip(1) * uvNext(end-1) + ip(2) * uvNext(end) + ip(3) * wvNext(1) + ip(4) * wvNext(2)); wvNext];
%             wv = [(ip(1) * uv(end-1) + ip(2) * uv(end) + ip(3) * wv(1) + ip(4) * wv(2)); wv];
            wvNext = [wvNextmh; wvNext];
            wv = [wvmh; wv];
            wvmh =  fliplr(customIp) * [uv(end-1:end); wv(1:2)];
            wpNext = [(customIp(4) * upNext(end-1) + customIp(3) * upNext(end) + customIp(2) * wpNext(1) + customIp(1) * wpNext(2)); wpNext];
            wp = [(customIp(4) * up(end-1) + customIp(3) * up(end) + customIp(2) * wp(1) + customIp(1) * wp(2)); wp];
               
        end
        [S, SHalf, SBar] = setTube(N+1);
        % insert matrix creation here
        
        potScalingU = ones(length(up),1);
        potScalingW = ones(length(wp),1);
        potScalingU(1) = 0.5;
        potScalingW(end) = 0.5;
        potScalingU(end) = 0.5;
        potScalingW(1) = 0.5;
    end   
    
    % remove point if N^n < N^{n-1}
    if N < NPrev
        if mod(N,2) == 0
            uvMph =  uv(end);
            uvNext = uvNext(1:end-1);
            uv = uv(1:end-1);
            upNext = upNext(1:end-1);
            up = up(1:end-1);

        else 
            wvmh = wv(1);
            wvNext = wvNext(2:end);
            wv = wv(2:end);
            wpNext = wpNext(2:end);
            wp = wp(2:end);

        end
        [S, SHalf, SBar] = setTube(N+1);
        
        potScalingU = ones(length(up),1);
        potScalingW = ones(length(wp),1);
        potScalingU(1) = 0.5;
        potScalingW(end) = 0.5;
        potScalingU(end) = 0.5
        potScalingW(1) = 0.5;

    end
    upRange = 2:length(up)-1;         % range without boundaries
    wpRange = 2:length(wp)-1;

    alfSave(n) = alf;
    A = [1, -ip(4); ...
         -ip(4), 1];
    v = [ip(1) * up(end-2) + ip(2) * up(end-1) + ip(3) * up(end); ...
         ip(3) * wp(1) + ip(2) * wp(2) + ip(1) * wp(3)];
    solut = A \ v;
    
    %% Calculate velocities
    uvNext = uv - lambda / (rho * c) * (up(2:end) - up(1:end-1));
    uvNextMph = uvMph - lambda / (rho * c) * (solut(2) - up(end));
    
    wvNext = wv - lambda / (rho * c) * (wp(2:end) - wp(1:end-1));
    wvNextmh = wvmh - lambda / (rho * c) * (wp(1) - solut(1));
    
    %% Calculate pressure
    upNext(upRange) = up(upRange) - rho * c * lambda ./ SBar(upRange) .* (SHalf(upRange) .* uvNext(upRange) - SHalf(upRange-1) .* uvNext(upRange-1));
    upNext(end) = up(end) - rho * c * lambda ./ SBar(length(up)) .* (SHalf(length(up)) .* uvNextMph - SHalf(length(up) - 1) .* uvNext(end));
    
    wpNext(wpRange) = wp(wpRange) - rho * c * lambda ./ SBar(wpRange + length(up) - 1) .* (SHalf(wpRange + length(up) - 1) .* wvNext(wpRange) - SHalf(wpRange + length(up) - 2) .* wvNext(wpRange-1));
    wpNext(1) = wp(1) - rho * c * lambda ./ SBar(length(up)) .* (SHalf(length(up)) .* wvNext(1) - SHalf(length(up) - 1) .* wvNextmh);

    %% Set output from output position
    out(n) = wp(end-5);
    
    %% Energies
    kinEnergyU(n) = rho / 2 * h * sum(SHalf(1:length(uv)) .* uvNext .* uv);
    potEnergyU(n) = 1/(2 * rho * c^2) * h * sum (SBar(1:length(up)) .* potScalingU .* up.^2);
    totEnergyU(n) = kinEnergyU(n) + potEnergyU(n);
    
    kinEnergyW(n) = rho / 2 * h * sum(SHalf(length(uv)+1:end) .* wvNext .* wv);
    potEnergyW(n) = 1/(2 * rho * c^2) * h * sum (SBar(length(up):end) .* potScalingW .* wp.^2);
    totEnergyW(n) = kinEnergyW(n) + potEnergyW(n);

    totEnergy(n) = totEnergyU(n) + totEnergyW(n);

    %% Draw things
    if drawThings && mod (n, drawSpeed) == 0 && n > drawStart
        hLocsLeft = (0:(length(up))-1) * h;
        hLocsRight = flip(L - ((0:(length(wp)-1)) * h));        
        
        % Plot the velocity
        subplot(4,1,1)
        cla
        hold on;
        plot(hLocsLeft / L, up, '-o');
        plot(hLocsRight / L, wp, '-o');
        plot(hLocsLeft(1:end) / L + 0.5 / N, [uvNext; uvNextMph] * 100, 'Marker', '.', 'MarkerSize', 10);
        plot(hLocsRight(1:end) / L - 0.5 / N, [wvNextmh; wvNext] * 100, 'Marker', '.', 'MarkerSize', 10);
%         xlim([0.4 0.6])
        ylim([-1 1])
        % Plot scaled energy
        subplot(4,1,2)
        hold off
        plot(kinEnergyU(1:n) + kinEnergyW(1:n))
        hold on;
        plot(potEnergyU(1:n) + potEnergyW(1:n))
%         plot(totEnergy(1:n) / totEnergy(1) - 1)
        subplot(4,1,3)
        plot(totEnergy(1:n) - totEnergy(1))
        
        subplot(4,1,4)
        plot(alfSave(1:n))
        pause(0.05);
        drawnow;
        
    end

    %% Update states
    uv = uvNext;
    up = upNext;
        
    wv = wvNext;
    wp = wpNext;

    uvMph = uvNextMph;
    wvmh = wvNextmh;

end   

function [S, SHalf, SBar] = setTube(N)
%     mp = linspace(0.001, 0.001, floor(N/40));           % mouthpiece
%     m2t = linspace(mp(end), 0.0005, floor(N/40));       % mouthpiece to tube
%     
%     alpha = 0.25;                                       % bell
%     b = m2t(end) * exp(alpha * (0:18));
%     pointsLeft = N - length([mp, m2t, b]);
%     tube = linspace(m2t(end), m2t(end), pointsLeft);    % tube

%     S = [mp, m2t, tube, b]';                            % True geometry
    S = ones(N, 1);
    % Calculate approximations to the geometry
    SHalf = (S(1:N-1) + S(2:N)) * 0.5;                  % mu_{x+}
    SBar = (SHalf(1:end-1) + SHalf(2:end)) * 0.5;
    SBar = [S(1); SBar; S(end)];                        % mu_{x-}S_{l+1/2}
    
end
