%{
    Webster's equation time varying
%}

clear all;
close all;

% drawing variables
drawThings = true;
drawSpeed = 10;
centered = false;
dynamic = true;
damping = true;

impulse = true;

fs = 44100;         % Sample rate (Hz)
k = 1/fs;           % Time step (s)
lengthSound = fs * 2; % Duration (s)

c = 343;            % Wave speed (m/s)
h = c * k;          % Grid spacing (m)
L = 1;

N = floor(0.95 * L/h);  % Number of points (-)
h = L/N;                 % Recalculate gridspacing from number of points
lambdaSq = (c * k / h)^2
altLambdaSq = 2 * k^2 * c^2 / h^2;

%% operators

% averaging
muXF = 0.5 * (sparse(2:N, 1:N-1, ones(1, N-1), N, N) + ...
        sparse(1:N, 1:N, ones(1, N), N, N));
muXB = 0.5 * (sparse(1:N-1, 2:N, ones(1, N-1), N, N) + ...
        sparse(1:N, 1:N, ones(1, N), N, N));
muXX = 0.25 * (sparse(1:N-1, 2:N, ones(1, N-1), N, N) + ...
        sparse(1:N, 1:N, 2 * ones(1, N), N, N) + ...
        sparse(2:N, 1:N-1, ones(1, N-1), N, N));


% Set cross-sectional geometry
[S, SHalf, SBar] = setTubeOverTime (N, lengthSound, dynamic);

SNph = 2 * SBar(N, :) - SHalf(end, :);
SOnemh = 2 * SBar(1, :) - SHalf(1, :);

if damping
    a1 = 1 / (2 * (0.8216)^2 * c);
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

% Set rgs
rg = 2:N-1;          % "clamped"
potEnergyRange = 2:N-1;    % energy rg
% potEnergyRange = 1:N-1; % rg for the potential energy

% Right below NSS Eq. (5.28) Bilbao explains that a primed inner product is used when the boundary condition is centered
scaling = ones(N,1);
if centered
    scaling([1 N]) = 0.5;
end

muTT = 0.25 * [ones(N-1,1), 2 * ones(N-1,1), ones(N-1,1)];
for n = 2:lengthSound - 1

    muTTSHalf = sum(muTT .* SHalf(:, n-1:n+1), 2);
    
    muTTSNph = sum(muTT(1,:) .* SNph(n-1:n+1));
    muTTSOnemh = sum(muTT(1,:) .* SOnemh(n-1:n+1));
%     sum(muTTSmh - SHalf(1:end-1, n))
    % calculate scheme
    uNext(rg) = ((sum(2 * muTT(rg, :) .* SBar(rg, n-1:n+1), 2) - lambdaSq * (muTTSHalf(rg) + muTTSHalf(rg-1))) .* u(rg) ...
        + lambdaSq * muTTSHalf(rg) .* u(rg+1)...
        + lambdaSq * muTTSHalf(rg-1) .* u(rg-1)...
        - 1/2 * (SBar(rg, n) + SBar(rg, n-1)) .* uPrev(rg)) ./ (1/2 * (SBar(rg, n+1) + SBar(rg, n)));
    if centered
%         uNext(1) = 2 * (1 - lambdaSq) * u(1) - uPrev(1) + lambdaSq * 2 * u(2) + 2 * h * lambdaSq * SOnemh / SBar(1) * in(n);
        uNext(N) = ((sum(2 * muTT(end, :) .* SBar(end, n-1:n+1), 2) - lambdaSq * (muTTSNph + muTTSHalf(end))) .* u(N)...
            + lambdaSq * (muTTSNph + muTTSHalf(end)) * u(N-1)...
            + (h * lambdaSq * muTTSNph * (a1/k - a2) - 1/2 * (SBar(N, n) + SBar(N, n-1))) * uPrev(N))...
            / (1/2 * (SBar(N, n+1) + SBar(N, n)) + h * lambdaSq * (muTTSNph * (a1 / k + a2)));
    else
%         uNext(1) = 2 * (1 - lambdaSq) * u(1) - uPrev(1) + lambdaSq * SHalf(1) / SBar(1) * u(2) + lambdaSq * SOnemh / SBar(1) * u(1);
        uNext(N) = ((sum(2 * muTT(end, :) .* SBar(end, n-1:n+1), 2) - lambdaSq * (muTTSNph + muTTSHalf(end))) .* u(N)...
            + lambdaSq * muTTSNph * u(N)...
            + lambdaSq * muTTSHalf(end) * u(N-1)...
            + (h * lambdaSq * muTTSNph * (a1/(2*k) - a2/2) - 1/2 * (SBar(N, n) + SBar(N, n-1))) * uPrev(N))...
            / (1/2 * (SBar(N, n+1) + SBar(N, n)) + h * lambdaSq * muTTSNph * (a1 / (2*k) + a2/2));
    end

    % set output from output position
    out(n) = uNext(outputPos);
    
%     % energies
%     kinEnergy(n) = h * 1/2 * sum(SBar .* scaling .* (1/k * (u - uPrev)).^2);
%     potEnergy(n) = -h * c^2 / 2 * 1/h^2 * sum(SHalf(potEnergyrg)...
%         .* (u(potEnergyrg+1) - u(potEnergyrg)) .* (uPrev(potEnergyrg+1) - uPrev(potEnergyrg)));
%     if centered
%         boundaryEnergy(n) = c^2 * SNph * a2 / 2 * (u(N)^2 + uPrev(N)^2);
%     else
%         boundaryEnergy(n) = c^2 * SNph * a2 / 4 * (u(N)^2 + uPrev(N)^2);
%     end
%     totEnergy(n) = kinEnergy(n) - potEnergy(n) + boundaryEnergy(n);
%     
%     %% Rate-of-Changes of energy
    rOCkinEnergy(n) = h / (8 * k^3) * sum((uNext - uPrev) .* ((SBar(:, n+1) - SBar(:, n-1)) .* (uNext - uPrev) + (SBar(:, n+1) + 2 * SBar(:, n) + SBar(:, n-1)) .* (uNext - 2 * u + uPrev))); % .* scaling
    rOCpotEnergy(n) = c^2 * h / (8 * k * h^2) * sum((uNext(potEnergyRange) - uPrev(potEnergyRange)) ...
        .* (1/4 * (S(potEnergyRange + 1, n+1) - S(potEnergyRange - 1, n+1)...
        + 2 * S(potEnergyRange + 1, n) - 2 * S(potEnergyRange - 1, n)...
        + S(potEnergyRange + 1, n-1) - S(potEnergyRange - 1, n-1))...
        .* (u(potEnergyRange + 1) - u(potEnergyRange - 1))...
        + (SBar(potEnergyRange, n+1) + 2 * SBar(potEnergyRange, n) + SBar(potEnergyRange, n-1))...
        .* (u(potEnergyRange + 1) - 2 * u(potEnergyRange) + u(potEnergyRange - 1))...
        ));
    
    if centered
        % left boundary
        rOCboundaryEnergy(n) = -c^2 * SOnemh * 1 / (2 * k) * (uNext(1) - uPrev(1)) * 1/h * (u(1) - u(2));
        % right boundary
        rOCboundaryEnergy(n) = rOCboundaryEnergy(n) + c^2 * SNph * (-2 * a1 * (1/(2*k) * (uNext(N) - uPrev(N)))^2 - a2 / (2*k) * (uNext(N)^2 - uPrev(N)^2) - 1 / (2 * h * k) * (uNext(N) - uPrev(N)) * (u(N) - u(N-1)));

    else
        % no boundaryEnergy at the left boundary as u(1) = u(2) -> ... * (u(1)-u(2)) = 0
%         rOCboundaryEnergy(n) = c^2 * muTTSNph * (-a1 * (1/(2*k) * (uNext(N) - uPrev(N)))^2 - a2 / (4*k) * (uNext(N)^2 - uPrev(N)^2));
        rOCboundaryEnergy(n) = c^2 * muTTSNph * 1 / (2*k) * (uNext(N) - uPrev(N)) * (-a1 /(2*k) * (uNext(N) - uPrev(N)) - a2 / 2 * (uNext(N) + uPrev(N)));

    end
    
    rOCtotEnergy(n) = rOCkinEnergy(n) - rOCpotEnergy(n);
%     
    % draw things
    if drawThings && mod (n, drawSpeed) == 0
%         plot(uNext)
        % plot state
        subplot(2,1,1)
        cla
        hold on;
        plotPressurePotential (uNext * 1/amp, S(:,n) * amp);
        plot(uNext)
        plot(S(:,n) * amp, 'k');
        plot(-S(:,n) * amp, 'k');
        xlim([1 N]);
        ylim([-max(S(:,n)) max(S(:,n))] * amp * 1.1);
        title("Pressure potential. n = " + num2str(n))
% 
%         % plot total energy
        subplot(2,1,2)
%         plot(rOCkinEnergy(1:n));
%         hold on;
%         plot(rOCpotEnergy(1:n));
%         hold off;
        plot(rOCtotEnergy(1:n))
        hold on;
        plot(rOCboundaryEnergy(1:n))
        hold off;
%         subplot(3,1,2)
%         if n > 10
%             hold off;
%             plot(totEnergy(10:n) / totEnergy(10) - 1);
% %             plot(totEnergy(10:n) - totEnergy(1));
% %             plot(kinBoundary(10:n))
% %             hold on;
% %             plot(potBoundEnergy(10:n))
%         end
%         title("Normalised energy (should be within machine precision)") 
%         
%         subplot(3,1,3)
%         plot(rOCtotEnergy(1:n))
% %         title("Rate of change of energy (should be 0-ish)") 
% %         hold off;
% %         plot(boundaryEnergyTest(1:n))
% %         hold on;
        drawnow;
%         
%         if makeVideo
%             M(frame) = getframe(gcf);
%             frame = frame + 1;
%         end
    end
   
    % update states
    uPrev = u;
    u = uNext;
end   
plot(out)

function [S, SHalf, SBar] = setTubeOverTime(N, lengthSound, dynamic)
    mp = ones(lengthSound, floor(N/20)) * 0.2;               % mouthpiece
    m2t = zeros(lengthSound, floor(N/20));
    b = zeros(lengthSound, 19);

    for i = 1:lengthSound
        m2t(i,:) = linspace(0.2, 0.1, floor(N/20));
    end
    alpha = 0.15;
    bPre = 0.1 * exp(alpha * (0:18));      % bell
    
    for i = 1:lengthSound
        b(i,:) = bPre;
    end
%     b = [linspace(0.1, 1, 18), 1];
    pointsLeft = N - length([mp(1,:), m2t(1,:), b(1,:)]);
    tube = b(1) * ones(lengthSound, pointsLeft);        % tube
    
    if dynamic
        bulgeWidth = floor(pointsLeft / 3);
        start = pointsLeft * (0.25 + 0.5 * (1:lengthSound)' / lengthSound) - bulgeWidth * 0.5;
        flooredStart = floor(start);
        raisedCos = 1 - (cos(2 * pi * ([1:bulgeWidth] - (start-floor(start))) / bulgeWidth) + 1) * 0.5;
        for i = 1:lengthSound
            tube(i, flooredStart(i):flooredStart(i)+bulgeWidth-1) = tube(i, flooredStart(i):flooredStart(i)+bulgeWidth-1) + raisedCos(i,:);
        end
    end
    S = [mp, m2t, tube, b]';
%     S = ones(N, lengthSound);
    SHalf = (S(1:N-1,:) + S(2:N,:)) * 0.5;            % mu_{x+}
    SBar = (SHalf(1:end-1,:) + SHalf(2:end,:)) * 0.5;
    SBar = [S(1) * ones(1, lengthSound); SBar; S(end) * ones(1, lengthSound)]; % mu_{x-}S_{l+1/2}
    
end
