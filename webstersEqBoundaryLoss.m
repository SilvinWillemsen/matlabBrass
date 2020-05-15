%{
    Webster's equation
%}

clear all;
close all;

% drawing variables
drawThings = true;
drawSpeed = 50;
centered = true;

impulse = false;

fs = 44100;         % Sample rate (Hz)
k = 1/fs;           % Time step (s)
lengthSound = fs*5; % Duration (s)

c = 343;            % Wave speed (m/s)
h = c * k;          % Grid spacing (m)
L = 1;              % Length

N = floor(L/h);             % Number of points (-)
h = L/N;                    % Recalculate gridspacing from number of points
lambdaSq = (c * k / h)^2    % Courant number

% Set cross-sectional geometry
[S, SHalf, SBar] = setTube (N);

a1 = 1 / (2 * (0.8216)^2 * c);              % loss term
a2 = L / (0.8216 * sqrt(S(1)*S(N)/pi));     % inertia coefficient
% a1 = 0;
% a2 = 0;

%Initialise states
uNext = zeros(N, 1);
u = zeros(N, 1);

amp = 1e-5;
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
range = 2:N-1;          % range without boundaries
potEnergyRange = 1:N-1; % range for the potential energy

% set up figure
figure(1);
subplot(2,1,1)
plot(sqrt(S), 'k');
hold on;
plot(-sqrt(S), 'k');
xlim([1 N]);
ylim([-max(sqrt(S)) max(sqrt(S))] * 1.1);
title("Pressure potential")

subplot(2,1,2)
title("Normalised energy (should be within machine precision)") 

% Problem 9.5 (weighted boundary conditions)
epsilonL = SHalf(1)/SBar(1);
epsilonR = SHalf(end)/SBar(N);

scaling = ones(N,1);
if centered
    scaling(1) = epsilonL / 2;
    scaling(N) = epsilonR / 2;
end

SNph = 2 * SBar(N) - SHalf(end);
SOnemh = 2 * SBar(1) - SHalf(1);

for n = 1:lengthSound
    
    % calculate scheme
    uNext(range) = 2 * (1 - lambdaSq) * u(range) - uPrev(range) + lambdaSq * ((SHalf(range) ./ SBar(range)) .* u(range+1) + (SHalf(range - 1) ./ SBar(range)) .* u(range-1));
    
    if centered
        uNext(1) = 2 * (1 - lambdaSq) * u(1) - uPrev(1) + lambdaSq * 2 * u(2) + 2 * h * lambdaSq * SOnemh / SBar(1) * in(n);
        uNext(N) = (2 * (1 - lambdaSq) * u(N) - uPrev(N) + lambdaSq * 2 * u(N-1) + h * lambdaSq * SNph / SBar(N) * (a1/k - a2) * uPrev(N)) / (1 + lambdaSq * SNph / SBar(N) * h * (a1/k + a2));
    else
        uNext(1) = 2 * (1 - lambdaSq) * u(1) - uPrev(1) + lambdaSq * SHalf(1) / SBar(1) * u(2) + lambdaSq * SOnemh / SBar(1) * u(1);
        uNext(N) = (2 * (1 - lambdaSq) * u(N) - uPrev(N) + lambdaSq * SNph * h / (2 * SBar(N)) * (a1 / k - a2) * uPrev(N) + lambdaSq * SNph / SBar(N) * u(N) + lambdaSq * SHalf(end) / SBar(N) * u(N-1)) / (1 + lambdaSq * SNph * h * (a1 / k + a2) / (2*SBar(N)));
    end

    % set output from output position
    out(n) = uNext(outputPos);
    
    % energies
    kinEnergy(n) = 1/2 * sum(h * SBar .* scaling .* (1/k * (u - uPrev)).^2);
    potEnergy(n) = -h * c^2 / 2 * 1/h^2 * sum(SHalf(potEnergyRange)...
        .* (u(potEnergyRange+1) - u(potEnergyRange)) .* (uPrev(potEnergyRange+1) - uPrev(potEnergyRange)));

    if centered
        boundaryEnergy(n) = 2 * (1 - epsilonR / 2) * SHalf(end) * c^2 * a2 / 4 * (u(N)^2 + uPrev(N)^2);
    else
        boundaryEnergy(n) = c^2 * SNph * a2 / 4 * (u(N)^2 + uPrev(N)^2);
    end
    
    totEnergy(n) = kinEnergy(n) - potEnergy(n) + boundaryEnergy(n);
    
    %% Rate-of-Changes of energy
    rOCkinEnergy(n) = h / (2*k^3) * sum(SBar .* scaling .* (uNext - 2 * u + uPrev) .* (uNext - uPrev));
    rOCpotEnergy(n) = -c^2 / (2 * k * h) * sum(SHalf .* (uNext(potEnergyRange+1) - uNext(potEnergyRange) - uPrev(potEnergyRange+1) + uPrev(potEnergyRange)) .* (u(potEnergyRange+1) - u(potEnergyRange)));
    
    if centered
        rOCboundaryEnergy(n) = -(2-epsilonL) * SHalf(1) * c^2 * 1/(2*k) * (uNext(1) - uPrev(1)) * (-in(n));
        rOCboundaryEnergy(n) = rOCboundaryEnergy(n) + (2-epsilonR) * SHalf(end) * c^2 * (-a1 * (1/(2*k) * (uNext(N) - uPrev(N)))^2 - a2/(4*k) * (uNext(N)^2 - uPrev(N)^2));
    else
        % no boundaryEnergy at the left boundary as u(1) = u(2) -> ... * (u(1)-u(2)) = 0
        rOCboundaryEnergy(n) = c^2 * SNph * (-a1 * (1/(2*k) * (uNext(N) - uPrev(N)))^2 - a2 / (4*k) * (uNext(N)^2 - uPrev(N)^2));
    end
    
    rOCtotEnergy(n) = rOCkinEnergy(n) - rOCpotEnergy(n) - rOCboundaryEnergy(n);
    
    % draw things
    if drawThings && mod (n, drawSpeed) == 0
        
        % plot state
        subplot(3,1,1)
        cla
        hold on;
        plotPressurePotential (uNext * 1/amp, sqrt(S) * amp);
        plot(uNext * 10000 * amp)
        plot(sqrt(S) * amp, 'k');
        plot(-sqrt(S) * amp, 'k');
        xlim([1 N]);
        ylim([-max(sqrt(S)) max(sqrt(S))] * amp * 1.1);
        title("Pressure potential. n = " + num2str(n))

        % plot total energy
        subplot(3,1,2)
        if n > 10
            if impulse && a1 == 0
                plot(totEnergy(10:n) / totEnergy(10) - 1);
                title("Normalised energy (should be within machine precision)") 
            else
                plot(totEnergy(10:n))
                title("Total energy") 
            end
        end
        
        subplot(3,1,3)
        plot(rOCtotEnergy(1:n))
        title("Rate-of-change of energy (should be very close to 0)");
        drawnow;
    end
   
    % update states
    uPrev = u;
    u = uNext;
end   
plot(out)

function [S, SHalf, SBar] = setTube(N)
    mp = linspace(0.0005, 0.0005, floor(N/20));       % mouthpiece
    m2t = linspace(0.0005, 0.0001, floor(N/20));     % mouthpiece to tube
    alpha = 0.25;
    b = m2t(end) * exp(alpha * (0:18));               % bell
    pointsLeft = N - length([mp, m2t, b]);
    tube = linspace(m2t(end), m2t(end), pointsLeft);        % tube

    S = [mp, m2t, tube, b]';                        % True geometry

    % Calculate approximations to the geometry
    SHalf = (S(1:N-1) + S(2:N)) * 0.5;           	% mu_{x+}
    SBar = (SHalf(1:end-1) + SHalf(2:end)) * 0.5;
    SBar = [S(1); SBar; S(end)];                    % mu_{x-}S_{l+1/2}

end
