%{
    First order brass
%}

clear all;
close all;

% drawing variables
drawThings = false;
drawSpeed = 10;
centered = true;

impulse = true;

fs = 44100;         % Sample rate (Hz)
k = 1/fs;           % Time step (s)
lengthSound = fs; % Duration (s)

T = 26.85;
[c, rho, eta, nu, gamma] = calcThermoDynConstants(T);

h = c * k;          % Grid spacing (m)
L = 1;              % Length

N = floor(L/h);             % Number of points (-)
h = L/N;                    % Recalculate gridspacing from number of points

lambda = c * k / h

a1 = 1; % loss
% Set cross-sectional geometry
[S, SHalf, SBar] = setTube (N);

%Initialise states
pNext = zeros(N, 1);
p = zeros(N, 1);
vNext = zeros(N-1, 1);
v = zeros(N-1, 1);

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
    p(floor(N / 3) - 5 : floor(N / 3) + 5) = amp*hann(11);
end

% output
out = zeros(lengthSound, 1);
outputPos = floor(1/5 * N);

% Set ranges
pRange = 2:N-1;          % range without boundaries
vRange = 1:N-1;

scaling = ones(N,1);
if centered
    scaling(1) = 1 / 2;
    scaling(N) = 1 / 2;
end

SNph = 2 * SBar(N) - SHalf(end);
SOnemh = 2 * SBar(1) - SHalf(1);

%% Initialise energies
kinEnergy = zeros(lengthSound, 1);
potEnergy = zeros(lengthSound, 1);
totEnergy = zeros(lengthSound, 1);

for n = 1:lengthSound
    %% Schemes
    
    % calculate schemes
    vNext(vRange) = v(vRange) - lambda / (rho * c) * (p(vRange+1) - p(vRange));
    pNext(pRange) = p(pRange) - rho * c * lambda ./ SBar(pRange) .* (SHalf(pRange) .* vNext(pRange) - SHalf(pRange-1) .* vNext(pRange-1));
    
%     vNext(end) = v(vRange) - lambda / (rho * c) * (p(vRange+1) - p(vRange));
    pNext(N) = p(N) - rho * c * lambda ./ SBar(N) .* (-2 * SHalf(end) .* vNext(end));
    
    % set output from output position
    out(n) = p(outputPos);
    
    %% Energies
    kinEnergy(n) = 1/(2 * rho * c^2) * h * sum(SBar .* scaling .* p.^2);
    potEnergy(n) = rho / 2 * h * sum(SHalf .* vNext .* v);
    totEnergy(n) = kinEnergy(n) + potEnergy(n);
    
    % draw things
    if drawThings && mod (n, drawSpeed) == 0
        subplot(3,1,1)
        cla
        hold on;
        plotPressurePotential (p * 1/amp, sqrt(S) * amp);
        plot(p)
        plot(sqrt(S) * amp, 'k');
        plot(-sqrt(S) * amp, 'k');
        xlim([1 N]);
        ylim([-max(sqrt(S)) max(sqrt(S))] * amp * 1.1);
        title("Pressure");
        
        subplot(3,1,2)
        plot(vNext);
        title("Particle Velocity")

        subplot(3,1,3)
        plot(totEnergy(1:n) / totEnergy(1) - 1);
        title("Normalised total energy (should be 0 within machine precision)")
        
        drawnow;
    end
   
    % update states
    v = vNext;
    p = pNext;
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
%     S = ones(N,1);
    % Calculate approximations to the geometry
    SHalf = (S(1:N-1) + S(2:N)) * 0.5;           	% mu_{x+}
    SBar = (SHalf(1:end-1) + SHalf(2:end)) * 0.5;
    SBar = [S(1); SBar; S(end)];                    % mu_{x-}S_{l+1/2}
    
end
