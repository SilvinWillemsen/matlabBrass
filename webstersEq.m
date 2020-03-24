%{
    Webster's equation
%}

clear all;
close all;

% drawing variables
makeVideo = false;
drawThings = true;
drawSpeed = 5;

impulse = true;

fs = 44100;         % Sample rate (Hz)
k = 1/fs;           % Time step (s)
lengthSound = fs*5; % Duration (s)

c = 343;            % Wave speed (m/s)
h = c * k;          % Grid spacing (m)
L = 1;
s0 = 0;
N = floor(L * 0.95/h);  % Number of points (-)
h = L/N;                 % Recalculate gridspacing from number of points
lambdaSq = (c * k / h)^2

% Set cross-sectional geometry
[S, SHalf, Sbar] = setTube (0, N, lengthSound);

if ~impulse  
    % input signal
    t = (0:lengthSound - 1) / fs;
    freq = 272;
    in = cos(2 * pi * freq * t) - 0.5;
    in = (in + abs(in)) / 2; % subplus

    rampLength = 1000; 
    env = [linspace(0, 1, rampLength), ones(1, lengthSound - rampLength)];
    in = in .* env;
    in = in - sum(in) / length(in);
end

%Initialise states
uNext = zeros(N, 1);
u = zeros(N, 1);
if impulse
    in = zeros(lengthSound, 1);
    u(floor(N / 2) - 5 : floor(N / 2) + 5) = hann(11);
end
uPrev = u;

% output
out = zeros(lengthSound, 1);
outputPos = floor(1/5 * N);

% Initialise energy vectors
kinEnergy = zeros(lengthSound, 1);
potEnergy = zeros(lengthSound, 1);
boundaryEnergy = zeros(lengthSound, 1);
totEnergy = zeros(lengthSound, 1);

rOCkinEnergy = zeros(lengthSound, 1);
rOCpotEnergy = zeros(lengthSound, 1);
rOCtotEnergy = zeros(lengthSound, 1);

% Set ranges
range = 2:N-1;          % "clamped"
energyRange = 2:N-1;    % energy range

% set up figure
figure(1);
subplot(2,1,1)
plot(S, 'k');
hold on;
plot(-S, 'k');
xlim([1 N]);
ylim([-max(S) max(S)] * 1.1);
title("Pressure potential")

subplot(2,1,2)
title("Normalised energy (should be within machine precision)") 

% figure('Renderer', 'painters', 'Position', [200 200 600 200])

if makeVideo
    loops = fs / drawSpeed + 1;
    M(100) = struct('cdata',[],'colormap',[]);
    set(gca, 'Color', 'white')
    frame = 1;
end

for n = 1:lengthSound
    [S, SHalf, Sbar] = setTube(n, N, lengthSound);
    
    % calculate scheme
    uNext(range) = (2 * u(range) - uPrev(range) + lambdaSq * ((SHalf(range) ./ Sbar(range-1)) .* u(range+1) + (SHalf(range - 1) ./ Sbar(range - 1)) .* u(range-1) - 2 * u(range))...
        + s0 * k ./ Sbar(range - 1) .* uPrev(range)) ./ (1 + s0 * k ./ Sbar(range-1));
    uNext(1) = (2 * lambdaSq * u(2) + 2 * (1 - lambdaSq) * u(1) - uPrev(1) + 2 * h * SHalf(1) / (Sbar(1)) * in(n) + s0 * k / Sbar(1) * uPrev(1)) / (1 + s0 * k / Sbar(1));
    % set output from output position
    out(n) = uNext(outputPos);
    
    % energies
    kinEnergy(n) = h * 1/2 * sum(Sbar(energyRange - 1) .* (1/k * (u(energyRange) - uPrev(energyRange))).^2);
    potEnergy(n) = h * c^2 / 2 * 1/h^2 * (sum(SHalf(energyRange) .* (u(energyRange+1) - u(energyRange)) .* (uPrev(energyRange+1) - uPrev(energyRange))));
    
    if range(1) == 2
        kinEnergy(n) = kinEnergy(n) + h * 1/4 * sum(Sbar(1) * (1/k * (u(1) - uPrev(1)))^2);
        potEnergy(n) = potEnergy(n) + h * c^2 / 2 * 1/h^2 * (SHalf(1) .* (u(2) - u(1)) .* (uPrev(2) - uPrev(1)));
    end
    
    totEnergy(n) = kinEnergy(n) + potEnergy(n) + boundaryEnergy(n);
    
    % draw things
    if drawThings && mod (n, drawSpeed) == 0
        
        % plot state
        subplot(2,1,1)
        cla;
        hold on;
        plotPressurePotential (uNext, S);
        plot(uNext)
        plot(S, 'k');
        hold on;
        plot(-S, 'k');
        xlim([1 N]);
        ylim([-max(S) max(S)] * 1.1);
%         ylim([-7, 7]);
        title("Pressure potential. n = " + num2str(n))

%         plot total energy
        subplot(2,1,2)
        if n > 10
            plot(totEnergy(10:n) / totEnergy(10) - 1);
        end
%         hold off;
%         plot(kinEnergy(1:n) + potEnergy(1:n))
%         hold on;
%         plot(2 * boundaryEnergy(1:n))
    title("Normalised energy (should be within machine precision)") 
        drawnow;
        
        if makeVideo
            M(frame) = getframe(gcf);
            frame = frame + 1;
        end
    end
   
    % update states
    uPrev = u;
    u = uNext;
end   

plot(out)

function [S, SHalf, Sbar] = setTube(n, N, lengthSound)
    mp = linspace(0.2, 0.2, floor(N/20));               % mouthpiece
    m2t = linspace(0.2, 0.1, floor(N/20));              % mouthpiece to tube
    tube = linspace(0.1, 0.1, floor(3 * N / 4));        % tube
    bulgeWidth = floor(length(tube) / 3);
    start = length(tube) * (0.25 + 0.5 * n / lengthSound) - bulgeWidth * 0.5;
    flooredStart = floor(start);
    raisedCos = 1 - (cos(2 * pi * ([1:bulgeWidth] - (start-floor(start))) / bulgeWidth) + 1) * 0.5;
    tube(flooredStart:flooredStart+bulgeWidth-1) = tube(flooredStart:flooredStart+bulgeWidth-1) + raisedCos;
    
    pointsLeft = N - length([mp, m2t, tube]);
    alpha = 0.5 * n / lengthSound - 0.25;
    b = tube(end) * exp(alpha * (0:pointsLeft-1));      % bell
    S = [mp, m2t, tube, b]';
    S = ones(N, 1) * 1;
    % Calculate approximations to the geometry
    SHalf = (S(1:N-1) + S(2:N)) * 0.5;            % mu_{x+}
    Sbar = (SHalf(1:end-1) + SHalf(2:end)) * 0.5; % mu_{x-}S_{l+1/2}
end
