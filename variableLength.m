%{
    Variable length
%}

clear all;
close all;

% drawing variables
drawThings = true;
drawSpeed = 1000;
drawStart = 1;
centered = true;

impulse = true;

fs = 44100;         % Sample rate (Hz)
k = 1/fs;           % Time step (s)
lengthSound = fs*5; % Duration (s)

Ninit = 60;           % edit how many points you want
h = 1/Ninit;
cInit = h/k;
c = cInit;            % Wave speed (m/s)
h = cInit * k;          % Grid spacing (m)
L = 1;              % Length

N = floor(L/h);             % Number of points (-)
lambdaSq = (c * k / h)^2    % Courant number

% Set cross-sectional geometry
[S, SHalf, SBar] = setTube (N);

a1 = 1 / (2 * (0.8216)^2 * c);              % loss term
a2 = L / (0.8216 * sqrt(S(1)*S(N)/pi));     % inertia coefficient
a1 = 0;
a2 = 0;

%Initialise states
uNext = zeros(ceil(N/2) + 1, 1);
u = zeros(ceil(N/2) + 1, 1);

wNext = zeros(floor(N/2) + 1, 1);
w = zeros(floor(N/2) + 1, 1);

amp = 1e-5;
if ~impulse  
    % input signal
    t = (0:lengthSound - 1) / fs;
    freq = c/1.5155; 
    in = cos(2 * pi * freq * t) - 0.5;
    in = (in + abs(in)) / 2; % subplus
    in = in - sum(in) / length(in);
    in = in * amp / 10;
    rampLength = 1000; 
    env = [linspace(0, 1, rampLength), ones(1, lengthSound - rampLength)];
    in = in .* env;
    in = in - sum(in) / length(in);
else
    in = zeros(lengthSound, 1);
    u(floor(N / 3) - 5 : floor(N / 3) + 5) = amp*hann(11);
end
uPrev = u;
wPrev = w;
% output
out = zeros(lengthSound, 1);
outputPos = floor(1/5 * N);

% Initialise energy vectors
kinEnergyU = zeros(lengthSound, 1);
potEnergyU = zeros(lengthSound, 1);
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

scalingU = ones(length(u),1);
scalingW = ones(length(w),1);
if centered
    scalingU(1) = epsilonL / 2;
    scalingW(end) = epsilonR / 2;
    scalingU(end) = 0.5;
    scalingW(1) = 0.5;
end

SNph = 2 * SBar(N) - SHalf(end);
SOnemh = 2 * SBar(1) - SHalf(1);

eu = ones(length(u), 1);
Dxxu = spdiags([[SHalf(1:length(u)-1)./SBar(2:length(u)); 0] -2*eu ([0; 2; SHalf(2:length(u)-1)./SBar(2:length(u)-1)])], -1:1, length(u),length(u));
ew = ones(length(w), 1);
Dxxw = spdiags([([SHalf(length(u)-1:end-1)./SBar(length(u):end-1); 2; 0]) -2*ew ([0; SHalf(length(u)-1:end)./SBar(length(u)-1:end-1)])], -1:1, length(w),length(w));

interpolatedPoints = [0; 0];
changeC = true;
interpol = "cubic";
for n = 1:lengthSound
    % change wave speed
    if changeC
        c = cInit * (1-0.21 * n/fs);%* sin(0.5 * pi * n/fs));
    else
        c = c;
    end
    cSave(n) = c;
    
    % save previous state for comparison later
    NPrev = N;

    % recalculate gridspacing, points lambda^2 and alpha from new wave speed
    h = c*k;
    Ninit = 1/h;
    N = floor(1/h);
    Nsave(n) = N;
    hSave(n) = h;
    
    lambdaSq = c^2 * k^2 / h^2;

    alf = (Ninit - N);
  
    % calculate interpolator
    if interpol == "cubic"
        ip = [alf * (alf - 1) * (alf - 2) / -6, ...
                (alf - 1) * (alf + 1) * (alf - 2) / 2, ...
                alf * (alf + 1) * (alf - 2) / -2, ...
                alf * (alf + 1) * (alf - 1) / 6];
    else
        ip = [0, (1-alf), alf, 0];
    end

    if abs(N - NPrev) > 1
        disp('too fast')
    end
    
    % add point if N^n > N^{n-1}
    if N > NPrev
        if n > 55110
            disp("wait");
        end
        customIp = [alf * (alf + 1) / -((alf + 2) * (alf + 3)); ...
                        2 * alf * (alf + 1) / ((alf + 1) * (alf + 2)); ...
                        2 * (alf + 1) / ((alf + 2) * (alf + 1)); ...
                        2 * alf / -((alf + 3) * (alf + 2))];
        if mod(N,2) == 1
            uNext = [uNext; (customIp(1) * uNext(end-1) + customIp(2) * uNext(end) + customIp(3) * wNext(1) + customIp(4) * wNext(2))];
            u = [u; (customIp(1) * u(end-1) + customIp(2) * u(end) + customIp(3) * w(1) + customIp(4) * w(2))];
            uPrev = [uPrev; (customIp(1) * uPrev(end-1) + customIp(2) * uPrev(end) + customIp(3) * wPrev(1) + customIp(4) * wPrev(2))];
        else 
            wNext = [(customIp(4) * uNext(end-1) + customIp(3) * uNext(end) + customIp(2) * wNext(1) + customIp(1) * wNext(2)); wNext];
            w = [(customIp(4) * u(end-1) + customIp(3) * u(end) + customIp(2) * w(1) + customIp(1) * w(2)); w];
            wPrev = [(customIp(4) * uPrev(end-1) + customIp(3) * uPrev(end) + customIp(2) * wPrev(1) + customIp(1) * wPrev(2)); wPrev];
        end
        [S, SHalf, SBar] = setTube(N);
        % insert matrix creation here
        
        eu = ones(length(u), 1);
        Dxxu = spdiags([[SHalf(1:length(u)-1)./SBar(2:length(u)); 0] -2*eu ([0; 2; SHalf(2:length(u)-1)./SBar(2:length(u)-1)])], -1:1, length(u),length(u));
        ew = ones(length(w), 1);
        Dxxw = spdiags([([SHalf(length(u)-1:end-1)./SBar(length(u):end-1); 2; 0]) -2*ew ([0; SHalf(length(u)-1:end)./SBar(length(u)-1:end-1)])], -1:1, length(w),length(w));


        scalingU = ones(length(u),1);
        scalingW = ones(length(w),1);
        scalingU(1) = epsilonL / 2;
        scalingW(end) = epsilonR / 2;
        scalingU(end) = 0.5;
        scalingW(1) = 0.5;
    end   
    
    % remove point if N^n < N^{n-1}
    if N < NPrev
        if mod(N,2) == 0
            uNext = uNext(1:end-1);
            u = u(1:end-1);
            uPrev = uPrev(1:end-1);
        else 
            wNext = wNext(2:end);
            w = w(2:end);
            wPrev = wPrev(2:end);
        end
        [S, SHalf, SBar] = setTube(N);
        % insert matrix creation here
        eu = ones(length(u), 1);
        Dxxu = spdiags([[SHalf(1:length(u)-1)./SBar(2:length(u)); 0] -2*eu ([0; 2; SHalf(2:length(u)-1)./SBar(2:length(u)-1)])], -1:1, length(u),length(u));
        ew = ones(length(w), 1);
        Dxxw = spdiags([([SHalf(length(u)-1:end-1)./SBar(length(u):end-1); 2; 0]) -2*ew ([0; SHalf(length(u)-1:end)./SBar(length(u)-1:end-1)])], -1:1, length(w),length(w));
        
        scalingU = ones(length(u),1);
        scalingW = ones(length(w),1);
        scalingU(1) = epsilonL / 2;
        scalingW(end) = epsilonR / 2;
        scalingU(end) = 0.5;
        scalingW(1) = 0.5;

    end
    % calculate scheme
    interpolatedPointsPrev = interpolatedPoints;
    interpolatedPoints = [1, -ip(4); -ip(4), 1] \ [ip(3) * w(1) + ip(2) * w(2) + ip(1) * w(3);
                                        ip(1) * u(end-2) + ip(2) * u(end-1) + ip(3) * u(end)];

    range = 2:length(u) - 1;                                
%     uNext(range) = 2 * (1 - lambdaSq) * u(range) - uPrev(range) + lambdaSq * ((SHalf(range) ./ SBar(range)) .* u(range+1) + (SHalf(range - 1) ./ SBar(range)) .* u(range-1));
%     
%     uNext(1) = 2 * (1 - lambdaSq) * u(1) - uPrev(1) + lambdaSq * 2 * u(2);
%     uNext(end) = 2 * (1 - lambdaSq) * u(end) - uPrev(end) + lambdaSq * 2 * u(end-1);
    uNext = 2 * u + lambdaSq * Dxxu * u - uPrev;
    uNext(end) = uNext(end) + SHalf(length(u)) ./ SBar(length(u)) * lambdaSq * interpolatedPoints(1);
    
    uNext(1) = uNext(1) + 2 * h * lambdaSq * SOnemh / SBar(1) * in(n);
% %     uNext(end) = (2 * (1 - lambdaSq) * u(end) - uPrev(end) + lambdaSq * 2 * u(end-1) + h * lambdaSq * SNph / SBar(N) * (a1/k - a2) * uPrev(N)) / (1 + lambdaSq * SNph / SBar(N) * h * (a1/k + a2));
% 
    pressure = 1/(2*k) * (uNext(1) - uPrev(1)); 
    wNext = 2 * w + lambdaSq * Dxxw * w - wPrev;
    wNext(1) = wNext(1) + SHalf(length(u)-1) ./ SBar(length(u)-1) * lambdaSq * interpolatedPoints(2);
    
%     uNext(1) = uNext(1) + 2 * h * lambdaSq * SOnemh / SBar(1) * in(n);
%     wNext(end) = (wNext(end) + h * lambdaSq * SNph / SBar(N) * (a1/k - a2) * wPrev(end)) / (1 + lambdaSq * SNph / SBar(N) * h * (a1/k + a2));


    % set output from output position
    out(n) = uNext(outputPos);
    
    % energies
    
    kinEnergyU(n) = 1/2 * sum(h * SBar(1:length(u)) .* scalingU .* (1/k * (u - uPrev)).^2);
    potEnergyU(n) = -h * c^2 / 2 * 1/h^2 * sum(SHalf(1:length(u)-1)...
        .* (u(2:length(u)) - u(1:length(u)-1)) .* (uPrev(2:length(u)) - uPrev(1:length(u)-1)));
    
    totEnergyU(n) = kinEnergyU(n) - potEnergyU(n);
    
%     kinEnergyU(n) = 1/2 * sum(h * SBar(1:length(u)) .* scalingU .* (1/k * (u - uPrev)).^2);
%     potEnergyU(n) = -h * c^2 / 2 * 1/h^2 * sum(SHalf(1:length(u)-1)...
%         .* (u(2:length(u)) - u(1:length(u)-1)) .* (uPrev(2:length(u)) - uPrev(1:length(u)-1)));
%     potEnergyU(n) = potEnergyU(n) - h * c^2 / 2 * 1/h^2 * SHalf(length(u)) * (0 - u(length(u))) * (0 - uPrev(length(u)));
%    
%     totEnergyU(n) = kinEnergyU(n) - potEnergyU(n);
%     
    
    kinEnergyW(n) = 1/2 * sum(h * SBar(length(u)-1:end) .* scalingW .* (1/k * (w - wPrev)).^2);
    potEnergyW(n) = -h * c^2 / 2 * 1/h^2 * sum(SHalf(length(u)-1:end)...
        .* (w(2:length(w)) - w(1:length(w)-1)) .* (wPrev(2:length(w)) - wPrev(1:length(w)-1)));
    totEnergyW(n) = kinEnergyW(n) - potEnergyW(n);

    if centered
        boundaryEnergy(n) = 2 * (1 - epsilonR / 2) * SHalf(end) * c^2 * a2 / 4 * (u(end)^2 + uPrev(end)^2);
    else
        boundaryEnergy(n) = c^2 * SNph * a2 / 4 * (u(N)^2 + uPrev(N)^2);
    end
    
    totEnergy(n) = totEnergyU(n) + totEnergyW(n);%+ boundaryEnergy(n);
%     
%     %% Rate-of-Changes of energy
%     rOCkinEnergy(n) = h / (2*k^3) * sum(SBar .* scaling .* (uNext - 2 * u + uPrev) .* (uNext - uPrev));
%     rOCpotEnergy(n) = -c^2 / (2 * k * h) * sum(SHalf .* (uNext(potEnergyRange+1) - uNext(potEnergyRange) - uPrev(potEnergyRange+1) + uPrev(potEnergyRange)) .* (u(potEnergyRange+1) - u(potEnergyRange)));
%     
%     if centered
%         rOCboundaryEnergy(n) = -(2-epsilonL) * SHalf(1) * c^2 * 1/(2*k) * (uNext(1) - uPrev(1)) * (-in(n));
%         rOCboundaryEnergy(n) = rOCboundaryEnergy(n) + (2-epsilonR) * SHalf(end) * c^2 * (-a1 * (1/(2*k) * (uNext(N) - uPrev(N)))^2 - a2/(4*k) * (uNext(N)^2 - uPrev(N)^2));
%     else
%         % no boundaryEnergy at the left boundary as u(1) = u(2) -> ... * (u(1)-u(2)) = 0
%         rOCboundaryEnergy(n) = c^2 * SNph * (-a1 * (1/(2*k) * (uNext(N) - uPrev(N)))^2 - a2 / (4*k) * (uNext(N)^2 - uPrev(N)^2));
%     end
%     
%     rOCtotEnergy(n) = rOCkinEnergy(n) - rOCpotEnergy(n) - rOCboundaryEnergy(n);
    
    % draw things
    if drawThings && mod (n, drawSpeed) == 0 && n > drawStart
         hLocsLeft = (0:(length(u))-1) * h;
        hLocsRight = flip(1 - ((0:(length(w)-1)) * h));

        % plot state
        subplot(2,1,1)
        cla
        hold on;
        plotPressurePotential (uNext * 1/amp, sqrt(S(1:length(u))) * amp);
        plotPressurePotential (wNext * 1/amp, sqrt(S(length(u)-1:end)) * amp, length(u)-1);

        plot(hLocsLeft * N, uNext * 10000 * amp, 'LineWidth' , 1, 'Marker', '.', 'MarkerSize', 10, 'Color', 'b')
        plot(hLocsRight * N, wNext * 10000 * amp, 'Linewidth', 1,  'Marker', 'o', 'MarkerSize', 5, 'Color', 'r')

        plot(sqrt(S) * amp, 'k');
        plot(-sqrt(S) * amp, 'k');
        xlim([0 N]);
        ylim([-max(sqrt(S)) max(sqrt(S))] * amp * 1.1);
        title("Pressure potential. n = " + num2str(n))

        % plot total energy
        subplot(2,1,2)
        if n > 10
            if impulse && a1 == 0
                plot(totEnergy(10:n) / totEnergy(10) - 1);
                title("Normalised energy (should be within machine precision)") 
            else
                plot(totEnergy(10:n))
                title("Total energy") 
            end
        end
        
%         subplot(3,1,3)
%         plot(rOCtotEnergy(1:n))
%         title("Rate-of-change of energy (should be very close to 0)");
        drawnow;
    end
   
    % update states
    uPrev = u;
    u = uNext;
    
    wPrev = w;
    w = wNext;
    
end   
plot(out)

subplot(2,1,1)
hold off
plot(out)

subplot(2,1,2)
outfft = fft(out);
semilogy([0:lengthSound-1]'*fs/lengthSound, abs(outfft), 'r');
xlim([0 3*c])

function [S, SHalf, SBar] = setTube(N)
    mp = linspace(0.0005, 0.0005, floor(N/20));       % mouthpiece
    m2t = linspace(0.0005, 0.0001, floor(N/20));     % mouthpiece to tube
    alpha = 0.25;
    b = m2t(end) * exp(alpha * (0:18));               % bell
    pointsLeft = N - length([mp, m2t, b]);
    tube = linspace(m2t(end), m2t(end), pointsLeft);        % tube

    S = [mp, m2t, tube, b]';                        % True geometry
%     S = 0.001 * (1:length(S))';
%     S = 0.001 * ones(length(S), 1);
%     S = 0.1 * rand(length(S), 1);
%     load Ssave.mat
%     S = Ssave;
%     S(N/2-6:N/2+6) = 0.1;
%     S = linspace(0.001, 0.005, N)';
    % Calculate approximations to the geometry
    SHalf = (S(1:N-1) + S(2:N)) * 0.5;           	% mu_{x+}
    SBar = (SHalf(1:end-1) + SHalf(2:end)) * 0.5;
    SBar = [S(1); SBar; S(end)];                    % mu_{x-}S_{l+1/2}

end
