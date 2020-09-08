%{
    First order brass divided into two parts for dynamic changes
%}

clear all;
close all;

% drawing variables
drawThings = true;
drawSpeed = 10000;
centered = true;

fs = 44100;             % Sample rate (Hz)
k = 1/fs;               % Time step (s)
lengthSound = fs * 2;   % Duration (s)

%% viscothermal effects
T = 26.85;
[c, rho, eta, nu, gamma] = calcThermoDynConstants (T);

%% Tube variables
h = c * k;              % Grid spacing (m)
L = 1;                  % Length

N = floor(L/h);         % Number of points (-)
h = L/N;                % Recalculate gridspacing from number of points

lambda = c * k / h      % courant number

%% Set cross-sectional geometry
[S, SHalf, SBar] = setTube (N+1);

% Quick note: N is the number of spaces between the points so the number of points is N+1

%% Initialise states
upNext = zeros(ceil(N/2) + 1, 1);
up = zeros(ceil(N/2) + 1, 1);
uvNext = zeros(ceil(N/2) + 1, 1); 
uv = zeros(ceil(N/2) + 1, 1);

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

kinScalingU = ones(length(up),1);
kinScalingU(1) = 1 / 2;
kinScalingU(end) = 1 / 2;

kinScalingW = ones(length(wp),1);
kinScalingW(1) = 1 / 2;
kinScalingW(end) = 1 / 2;

potScalingU = ones(length(uv),1);
potScalingU(end) = 0.5;

potScalingW = ones(length(wv),1);
potScalingW(1) = 0.5;

for n = 1:lengthSound
    %% Calculate velocities before lip model
    uvNext(1:end-1) = uv(1:end-1) - lambda / (rho * c) * (up(2:end) - up(1:end-1));
    uvNext(end) = uv(end) - lambda / (rho * c) * (wp(2) - up(end));
    
    wvNext(2:end) = wv(2:end) - lambda / (rho * c) * (wp(3:end) - wp(2:end-1));
    wvNext(1) = wv(1) - lambda / (rho * c) * (wp(2) - wp(1));

    %% Calculate pressure
    upNext(upRange) = up(upRange) - rho * c * lambda ./ SBar(upRange) .* (SHalf(upRange) .* uvNext(upRange) - SHalf(upRange-1) .* uvNext(upRange-1));
    upNext(end) = up(end) - rho * c * lambda ./ SBar(length(up)) .* (SHalf(length(up)) .* uvNext(end) - SHalf(length(up) - 1) .* uvNext(length(up) - 1));
    
    wpNext(wpRange) = wp(wpRange) - rho * c * lambda ./ SBar(wpRange + length(up) - 1) .* (SHalf(wpRange + length(up) - 1) .* wvNext(wpRange) - SHalf(wpRange + length(up) - 2) .* wvNext(wpRange-1));
    wpNext(1) = wp(1) - rho * c * lambda ./ SBar(length(up)) .* (SHalf(length(up)) .* wvNext(1) - SHalf(length(up) - 1) .* uvNext(length(up) - 1));

    %% Set output from output position
    out(n) = wp(end-5);
    
    %% Energies
    kinEnergyU(n) = 1/(2 * rho * c^2) * h * sum (SBar(1:length(up)) .* kinScalingU .* up.^2);
    potEnergyU(n) = rho / 2 * h * sum(SHalf(1:length(uv)) .* potScalingU .* uvNext .* uv);
    totEnergyU(n) = kinEnergyU(n) + potEnergyU(n);
    
    kinEnergyW(n) = 1/(2 * rho * c^2) * h * sum (SBar(length(up):end) .* kinScalingW .* wp.^2);
    potEnergyW(n) = rho / 2 * h * sum(SHalf(length(uv):end) .* potScalingW .* wvNext .* wv);
    totEnergyW(n) = kinEnergyW(n) + potEnergyW(n);

    totEnergy(n) = totEnergyU(n) + totEnergyW(n);

    %% Draw things
    if drawThings && mod (n, drawSpeed) == 0
        hLocsLeft = (0:(length(up))-1) * h;
        hLocsRight = flip(L - ((0:(length(wp)-1)) * h));        
        
        % Plot the velocity
        subplot(2,1,1)
        cla
        hold on;
        plot(hLocsLeft * N / L, up, '-o');
        plot(hLocsRight * N / L, wp, '-o');
        plot(hLocsLeft * N / L + 0.5, uvNext * 100, 'Marker', '.', 'MarkerSize', 10);
        plot(hLocsRight(2:end) * N / L -0.5, wvNext * 100, 'Marker', '.', 'MarkerSize', 10);

        % Plot scaled energy
        subplot(2,1,2)
        plot(totEnergy(1:n) / totEnergy(1) - 1)
        pause(0.1);
        drawnow;
        
    end

    %% Update states
    uv = uvNext;
    up = upNext;
        
    wv = wvNext;
    wp = wpNext;

end   

function [S, SHalf, SBar] = setTube(N)
    mp = linspace(0.001, 0.001, floor(N/40));           % mouthpiece
    m2t = linspace(mp(end), 0.0005, floor(N/40));       % mouthpiece to tube
    
    alpha = 0.25;                                       % bell
    b = m2t(end) * exp(alpha * (0:18));
    pointsLeft = N - length([mp, m2t, b]);
    tube = linspace(m2t(end), m2t(end), pointsLeft);    % tube

    S = [mp, m2t, tube, b]';                            % True geometry

    % Calculate approximations to the geometry
    SHalf = (S(1:N-1) + S(2:N)) * 0.5;                  % mu_{x+}
    SBar = (SHalf(1:end-1) + SHalf(2:end)) * 0.5;
    SBar = [S(1); SBar; S(end)];                        % mu_{x-}S_{l+1/2}
    
end
