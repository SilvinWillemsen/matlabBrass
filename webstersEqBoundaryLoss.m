%{
    Webster's equation
%}

clear all;
close all;

% drawing variables
makeVideo = false;
drawThings = true;
drawSpeed = 6;
centered = true;
dynamic = false;
damping = true;

impulse = true;

fs = 44100;         % Sample rate (Hz)
k = 1/fs;           % Time step (s)
lengthSound = fs*1; % Duration (s)

c = 343;            % Wave speed (m/s)
h = c * k;          % Grid spacing (m)
L = 1;

N = floor(0.95 * L/h);  % Number of points (-)
h = L/N;                 % Recalculate gridspacing from number of points
lambdaSq = (c * k / h)^2

% Set cross-sectional geometry
[S, SHalf, SBar] = setTube (0, N, lengthSound, dynamic);

if damping
    a1 = L / (2 * (0.8216)^2 * c);
    a2 = L / (0.8216 * sqrt(S(1)*S(N)/pi));
    a1 = 0;
%     a2 = 0;
else
    a1 = 0;
    a2 = 0;
end

%Initialise states
uNext = zeros(N, 1);
u = zeros(N, 1);

amp = 1e-3;
if ~impulse  
    % input signal
    t = (0:lengthSound - 1) / fs;
    freq = 446/4;
    in = cos(2 * pi * freq * t) - 0.5;
    in = (in + abs(in)) / 2; % subplus
    in = in - sum(in) / length(in);
    in = in * amp;
    rampLength = 1000; 
    env = [linspace(0, 1, rampLength), ones(1, lengthSound - rampLength)];
    in = in .* env;
    in = in - sum(in) / length(in);
else
    in = zeros(lengthSound, 1);
    u(floor(N / 3) - 5 : floor(N / 3) + 5) = amp*hann(11);
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
rOCboundaryEnergy = zeros(lengthSound, 1);

% Set ranges
range = 2:N-1;          % "clamped"
energyRange = 2:N-1;    % energy range
potEnergyRange = 1:N-1; % range for the potential energy

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

if makeVideo
    loops = fs / drawSpeed + 1;
    M(100) = struct('cdata',[],'colormap',[]);
    set(gca, 'Color', 'white')
    frame = 1;
end

% Right below NSS Eq. (5.28) Bilbao explains that a primed inner product is used when the boundary condition is centered
scaling = ones(N,1);
if centered
    scaling([1 N]) = 0.5;
end
SNph = 2 * SBar(N) - SHalf(end);
SNmh = 2 * SBar(1) - SHalf(1);

for n = 1:lengthSound
    if dynamic
        [S, SHalf, SBar] = setTube(n, N, lengthSound, dynamic);
        SNph = 2 * SBar(N) - SHalf(end);
        SNmh = 2 * SBar(1) - SHalf(1);
    end
    
    % calculate scheme
    uNext(range) = 2 * (1 - lambdaSq) * u(range) - uPrev(range) + lambdaSq * ((SHalf(range) ./ SBar(range)) .* u(range+1) + (SHalf(range - 1) ./ SBar(range)) .* u(range-1));
    
    if centered
        uNext(1) = 2 * (1 - lambdaSq) * u(1) - uPrev(1) + lambdaSq * 2 * u(2) + 2 * h * lambdaSq * SNmh / SBar(1) * in(n);
        uNext(N) = (2 * (1 - lambdaSq) * u(N) - uPrev(N) + lambdaSq * 2 * u(N-1) + h * lambdaSq * SNph / SBar(N) * (a1/k - a2) * uPrev(N)) / (1 + lambdaSq * SNph / SBar(N) * h * (a1/k + a2));
    else
        uNext(1) = 2 * (1 - lambdaSq) * u(1) - uPrev(1) + lambdaSq * SHalf(1) / SBar(1) * u(2) + lambdaSq * SNmh / SBar(1) * u(1);
        uNext(N) = (2 * (1 - lambdaSq) * u(N) - uPrev(N) + lambdaSq * SNph * h / (2 * SBar(N)) * (a1 / k - a2) * uPrev(N) + lambdaSq * SNph / SBar(N) * u(N) + lambdaSq * SHalf(end) / SBar(N) * u(N-1)) / (1 + lambdaSq * SNph * h * (a1 / k + a2) / (2*SBar(N)));
    end

    % set output from output position
    out(n) = uNext(outputPos);
    
    % energies
    kinEnergy(n) = h * 1/2 * sum(SBar .* scaling .* (1/k * (u - uPrev)).^2);
    potEnergy(n) = -h * c^2 / 2 * 1/h^2 * sum(SHalf(potEnergyRange)...
        .* (u(potEnergyRange+1) - u(potEnergyRange)) .* (uPrev(potEnergyRange+1) - uPrev(potEnergyRange)));
    if centered
        boundaryEnergy(n) = c^2 * SNph * a2 / 2 * (u(N)^2 + uPrev(N)^2);
    else
        boundaryEnergy(n) = c^2 * SNph * a2 / 4 * (u(N)^2 + uPrev(N)^2);
    end
    totEnergy(n) = kinEnergy(n) - potEnergy(n) + boundaryEnergy(n);
    
    %% Rate-of-Changes of energy
    rOCkinEnergy(n) = h / (2 * k^3) * sum(SBar .* (uNext - 2 * u + uPrev) .* (uNext - uPrev)); % .* scaling
    rOCpotEnergy(n) = -c^2 / (2 * k * h) * sum(SHalf .* (uNext(potEnergyRange+1) - uNext(potEnergyRange) - uPrev(potEnergyRange+1) + uPrev(potEnergyRange)) .* (u(potEnergyRange+1) - u(potEnergyRange)));
    
    if centered
        % left boundary
        rOCboundaryEnergy(n) = -c^2 * SNmh * 1 / (2 * k) * (uNext(1) - uPrev(1)) * 1/h * (u(1) - u(2));
        % right boundary
        rOCboundaryEnergy(n) = rOCboundaryEnergy(n) + c^2 * SNph * (-2 * a1 * (1/(2*k) * (uNext(N) - uPrev(N)))^2 - a2 / (2*k) * (uNext(N)^2 - uPrev(N)^2) - 1 / (2 * h * k) * (uNext(N) - uPrev(N)) * (u(N) - u(N-1)));

    else
        % no boundaryEnergy at the left boundary as u(1) = u(2) -> ... * (u(1)-u(2)) = 0
% correct: 
        rOCboundaryEnergy(n) = c^2 * SNph * (-a1 * (1/(2*k) * (uNext(N) - uPrev(N)))^2 - a2 / (4*k) * (uNext(N)^2 - uPrev(N)^2));
% % stefan's book
%         boundaryEnergy(n) = c^2 * SNph * (-a1 * (1/(2*k) * (uNext(N) - uPrev(N)))^2 - a2 / (8*k) * (uNext(N)^2 + 2 * uNext(N) * u(n) - 2 * u(N) * uPrev(N) - uPrev(N)^2));
    end
    
    rOCtotEnergy(n) = rOCkinEnergy(n) - rOCpotEnergy(n) - rOCboundaryEnergy(n);
    
    % draw things
    if drawThings && mod (n, drawSpeed) == 0
        
        % plot state
        subplot(3,1,1)
        cla
        hold on;
        plotPressurePotential (uNext * 1/amp, S * amp);
        plot(uNext)
        plot(S * amp, 'k');
        plot(-S * amp, 'k');
        xlim([1 N]);
        ylim([-max(S) max(S)] * amp * 1.1);
        title("Pressure potential. n = " + num2str(n))

        % plot total energy
        subplot(3,1,2)
        if n > 10
            hold off;
            plot(totEnergy(10:n) / totEnergy(10) - 1);
%             plot(totEnergy(10:n) - totEnergy(1));
%             plot(kinBoundary(10:n))
%             hold on;
%             plot(potBoundEnergy(10:n))
        end
        title("Normalised energy (should be within machine precision)") 
        
        subplot(3,1,3)
        plot(rOCtotEnergy(1:n))
%         title("Rate of change of energy (should be 0-ish)") 
%         hold off;
%         plot(boundaryEnergyTest(1:n))
%         hold on;
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
%     b = [linspace(0.1, 1, 18), 1];
    pointsLeft = N - length([mp, m2t, b]);
    tube = linspace(b(1), b(1), pointsLeft);        % tube
    
    if dynamic
        bulgeWidth = floor(length(tube) / 3);
        start = length(tube) * (0.25 + 0.5 * n / lengthSound) - bulgeWidth * 0.5;
        flooredStart = floor(start);
        raisedCos = 1 - (cos(2 * pi * ([1:bulgeWidth] - (start-floor(start))) / bulgeWidth) + 1) * 0.5;
        tube(flooredStart:flooredStart+bulgeWidth-1) = tube(flooredStart:flooredStart+bulgeWidth-1) + raisedCos;
    end
    S = [mp, m2t, tube, b]';
%     S = [ones(floor(N/2),1); 2 * ones(ceil(N/2),1)];
%     S(1) = 1.1;
%     S(N) = 4;
%     S = ones(N, 1) * 1;
    % Calculate approximations to the geometry
%     S = 0.1 + rand(N, 1) * 0.03;
%     S = linspace(1, 2, N)';
    SHalf = (S(1:N-1) + S(2:N)) * 0.5;            % mu_{x+}
    SBar = (SHalf(1:end-1) + SHalf(2:end)) * 0.5;
    SBar = [S(1); SBar; S(end)]; % mu_{x-}S_{l+1/2}
    
%     SBar = [SBar(1); SBar; SBar(end)]; % mu_{x-}S_{l+1/2}

end
