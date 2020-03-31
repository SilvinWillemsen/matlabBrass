%{
    Webster's equation
%}

clear all;
close all;

% drawing variables
makeVideo = false;
drawThings = true;
drawSpeed = 5;
centered = false;
dynamic = false;

impulse = true;

fs = 44100;         % Sample rate (Hz)
k = 1/fs;           % Time step (s)
lengthSound = fs*1; % Duration (s)

c = 200;            % Wave speed (m/s)
h = c * k;          % Grid spacing (m)
L = 1;

N = floor(0.95 * L/h);  % Number of points (-)
h = L/N;                 % Recalculate gridspacing from number of points
lambdaSq = (c * k / h)^2

% Set cross-sectional geometry
[S, SHalf, SBar] = setTube (0, N, lengthSound, dynamic);

a1 = L / (2 * (0.8216)^2 * c);
a2 = L / (0.8216 * sqrt(S(1)*S(N)/pi));
a1 = 0;
a2 = 0;
if ~impulse  
    % input signal
    t = (0:lengthSound - 1) / fs;
    freq = 446/4;
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
    u(floor(N * 2 / 3) - 5 : floor(N * 2 / 3) + 5) = 1e-3*hann(11);
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
scaling = ones(N,1);
scaling([1 N]) = 2;

for n = 1:lengthSound
    [S, SHalf, SBar] = setTube(n, N, lengthSound, dynamic);
    
    % calculate scheme
    uNext(range) = 2 * u(range) - uPrev(range) + lambdaSq * ((SHalf(range) ./ SBar(range-1)) .* u(range+1) + (SHalf(range - 1) ./ SBar(range - 1)) .* u(range-1) - 2 * u(range));
    if centered
        uNext(1) = 2 * lambdaSq * u(2) + 2 * (1 - lambdaSq) * u(1) - uPrev(1) + 2 * h * lambdaSq * SHalf(1) / SBar(1) * in(n);
        uNext(N) = (2 * lambdaSq * u(N-1) + 2 * (1 - lambdaSq) * u(N) - uPrev(N) + 2 * h * lambdaSq * (a1/k - a2) * uPrev(N)) / (1 + lambdaSq * 2 * h * (a1/k + a2));
    else
        uNext(1) = 2 * (1 - lambdaSq) * u(1) - uPrev(1) + 2 * h * lambdaSq * (SHalf(1) / SBar(1)) * in(n) + lambdaSq * SHalf(2) / SBar(1) * u(2) + lambdaSq * (SHalf(1) / SBar(1)) * u(1);
        uNext(N) = (2 * (1 - lambdaSq) * u(N) - uPrev(N) + 2 * h * lambdaSq * (a1/k - a2) * uPrev(N) + lambdaSq * SHalf(end-1) / SBar(end) * u(N-1) + lambdaSq * SHalf(end-1) / SBar(end) * u(N)) / (1 + lambdaSq * 2 * h * (a1/k + a2));
    end
    % set output from output position
    out(n) = uNext(outputPos);
    
    % energies
    kinEnergy(n) = h * 1/2 * sum(SBar(energyRange - 1) .* (1/k * (u(energyRange) - uPrev(energyRange))).^2);
    potEnergy(n) = h * c^2 / 2 * 1/h^2 * (sum(SHalf(energyRange) .* (u(energyRange+1) - u(energyRange)) .* (uPrev(energyRange+1) - uPrev(energyRange))));
    
    if range(1) == 2
        if centered
%             % SHalf works when the boundary is changing (i.e. \delta_xS_N
%             % != 0), but I don't know why..
            kinEnergy(n) = kinEnergy(n) + h * 1/4 * sum(SHalf(1) * (1/k * (u(1) - uPrev(1)))^2);
            kinEnergy(n) = kinEnergy(n) + h * 1/4 * sum(SHalf(end) * (1/k * (u(N) - uPrev(N)))^2);
        else
            kinEnergy(n) = kinEnergy(n) + h * 1/2 * sum(SHalf(1) * (1/k * (u(1) - uPrev(1)))^2);
            kinEnergy(n) = kinEnergy(n) + h * 1/2 * sum(SHalf(end) * (1/k * (u(N) - uPrev(N)))^2);
        end
            potEnergy(n) = potEnergy(n) + h * c^2 / 2 * 1/h^2 * (SHalf(1) .* (u(2) - u(1)) .* (uPrev(2) - uPrev(1)));
    end
    
    totEnergy(n) = kinEnergy(n) + potEnergy(n);
    
    
    %% Rate-of-Changes of energy
    rOCkinEnergy(n) = h / (2 * k^3) * sum(S .* (uNext - 2 * u + uPrev) .* (uNext - uPrev)); % .* scaling
    rOCpotEnergy(n) = h * c^2 / (2*k*h^2) * sum(SBar .* (u(energyRange+1) - 2 * u(energyRange) + u(energyRange-1)).* (uNext(energyRange) - uPrev(energyRange)));
    
    if centered
        rOCpotEnergy(n) = rOCpotEnergy(n) + h * c^2 / (2 * k * h^2) * SBar(1) * (2 * u(2) - 2 * u(1)) * (uNext(1) - uPrev(1));
        rOCpotEnergy(n) = rOCpotEnergy(n) + h * c^2 / (2 * k * h^2) * SBar(end) * (2 * u(N-1) - 2 * u(N)) * (uNext(N) - uPrev(N)); % to give the same scaling to the other boundary (?)
    else
        rOCpotEnergy(n) = rOCpotEnergy(n) + h * c^2 / (4 * k * h^2) * SBar(1) * (2 * u(2) - 2 * u(1)) * (uNext(1) - uPrev(1));
        rOCpotEnergy(n) = rOCpotEnergy(n) + h * c^2 / (4 * k * h^2) * SBar(end) * (2 * u(N-1) - 2 * u(N)) * (uNext(N) - uPrev(N)); % to give the same scaling to the other boundary (?)
    end
    
    boundaryEnergy(n) = c^2 * SHalf(end) * (-a1 * (1/k * (uNext(N) - uPrev(N)))^2 - a2 / (4*k) * (uNext(N)^2 - uPrev(N)^2));
%     boundaryEnergy2(n) = c^2 * SHalf(end) * (-a1 * (1/k * (uNext(N) - uPrev(N)))^2 - a2 / (4*k) * (uNext(N)^2 + 2 * uNext(N) * u(N) + 2 * u(N)*uPrev(N) - uPrev(N)^2));
    rOCtotEnergy(n) = rOCkinEnergy(n) - rOCpotEnergy(n);
%     boundaryEnergy1(n)-boundaryEnergy2(n)
    % draw things
    if drawThings && mod (n, drawSpeed) == 0
        
        % plot state
        subplot(3,1,1)
%         cla;
%         hold on;
hold off;
%         plotPressurePotential (uNext, S);
        plot(uNext)
%         plot(S, 'k');
%         hold on;
%         plot(-S, 'k');
%         xlim([1 N]);
%         ylim([-max(S) max(S)] * 1.1);
%         ylim([-7, 7]);
        title("Pressure potential. n = " + num2str(n))

%         plot total energy
        subplot(3,1,2)
        if n > 10
            plot(totEnergy(10:n) / totEnergy(10) - 1);
        end
        title("Normalised energy (should be within machine precision)") 
        
        subplot(3,1,3)
%         hold off;
%         plot(rOCkinEnergy(1:n))
%         hold on;
%         plot(rOCpotEnergy(1:n))
%         hold off;
        plot(rOCtotEnergy(1:n))
%         hold on;
%         plot(boundaryEnergy(1:n))
        title("Rate of change of energy (should be 0-ish)") 

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

function [S, SHalf, SBar] = setTube(n, N, lengthSound, dynamic)
    mp = linspace(0.2, 0.2, floor(N/20));               % mouthpiece
    m2t = linspace(0.2, 0.1, floor(N/20));              % mouthpiece to tube
    alpha = 0.15;
    b = 0.1 * exp(alpha * (0:18));      % bell
    pointsLeft = N - length([mp, m2t, b]);
    tube = linspace(b(1), b(1), pointsLeft);        % tube
    if dynamic
        bulgeWidth = floor(length(tube) / 3);
        start = length(tube) * (0.25 + 0.5 * n / lengthSound) - bulgeWidth * 0.5;
        flooredStart = floor(start);
        raisedCos = 1 - (cos(2 * pi * ([1:bulgeWidth] - (start-floor(start))) / bulgeWidth) + 1) * 0.5;
        tube(flooredStart:flooredStart+bulgeWidth-1) = tube(flooredStart:flooredStart+bulgeWidth-1) + raisedCos;
    end

%     S = [mp, m2t, tube, b]';
    S = ones(N, 1) * 1;
    % Calculate approximations to the geometry
    SHalf = (S(1:N-1) + S(2:N)) * 0.5;            % mu_{x+}
    SBar = (SHalf(1:end-1) + SHalf(2:end)) * 0.5; % mu_{x-}S_{l+1/2}
end
