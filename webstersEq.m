%{
    Webster's equation
%}

clear all;
close all;

% drawing variables
drawThings = true;
drawSpeed = 2;

fs = 44100;         % Sample rate (Hz)
k = 1/fs;           % Time step (s)
lengthSound = fs;   % Duration (s)

c = 343;            % Wave speed (m/s)
h = c * k;          % Grid spacing (m)
N = floor(0.95/h);  % Number of points (-)
h = 1/N;            % Recalculate gridspacing from number of points
lambdaSq = (c * k / h)^2

% Set cross-sectional geometry
S = ones(N,1);
SRange = floor(N/3):floor(2*N/3);
S(SRange) = S(SRange) - hann(length(SRange)) * 0.9;

% Calculate approximations to the geometry
SHalf = (S(1:N-1) + S(2:N)) * 0.5;            % mu_{x+}
Sbar = (SHalf(1:end-1) + SHalf(2:end)) * 0.5; % mu_{x-}S_{l+1/2}

%Initialise states
uNext = zeros(N, 1);
u = zeros(N, 1);

% .. with a hanning window
width = floor(N / 10);
loc = floor(N / 5);
u((1:width) + loc) = 2*hann(width);
uPrev = u;

% output
out = zeros(lengthSound, 1);
outputPos = floor(1/5 * N);

% Initialise energy vectors
kinEnergy = zeros(lengthSound, 1);
potEnergy = zeros(lengthSound, 1);
totEnergy = zeros(lengthSound, 1);

rOCkinEnergy = zeros(lengthSound, 1);
rOCpotEnergy = zeros(lengthSound, 1);
rOCtotEnergy = zeros(lengthSound, 1);

% Set ranges
range = 3:N-2;          % "clamped"
energyRange = 2:N-1;    % energy range

% set up figure
figure(1);
subplot(2,1,1)
plot(S, 'k');
hold on;
plot(-S, 'k');
xlim([1 N]);
ylim([-1.1 1.1]);
title("Pressure potential")
subplot(2,1,2)
title("Normalised energy (should be within machine precision)") 

for n = 1:lengthSound
    
    % calculate scheme
    uNext(range) = 2 * u(range) - uPrev(range) + lambdaSq * ((SHalf(range) ./ Sbar(range-1)) .* u(range+1) + (SHalf(range - 1) ./ Sbar(range - 1)) .* u(range-1) - 2 * u(range));
    
    % set output from output position
    out(n) = uNext(outputPos);
    
    % energies
    kinEnergy(n) = h * 1/2 * sum(Sbar(energyRange - 1) .* (1/k * (u(energyRange) - uPrev(energyRange))).^2);
    potEnergy(n) = h * c^2 / 2 * 1/h^2 * (sum(SHalf(energyRange) .* (u(energyRange+1) - u(energyRange)) .* (uPrev(energyRange+1) - uPrev(energyRange))));
    totEnergy(n) = kinEnergy(n) + potEnergy(n);
    
    % draw things
    if drawThings && mod (n, drawSpeed) == 0
        
        % plot state
        subplot(2,1,1)
        plotPressurePotential (uNext, S);
        
        % plot total energy
        subplot(2,1,2)
        if n > 1
            plot(totEnergy(1:n) / totEnergy(1) - 1);
        end
        title("Normalised energy (should be within machine precision)") 
        drawnow;
    end
   
    % update states
    uPrev = u;
    u = uNext;
end   

plot(out)