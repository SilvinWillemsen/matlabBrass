%{
    First order brass
%}

clear all;
close all;

% drawing variables
drawThings = false;
drawSpeed = 5;
centered = true;

impulse = true;

fs = 44100;         % Sample rate (Hz)
k = 1/fs;           % Time step (s)
lengthSound = fs * 20; % Duration (s)

%% Tube variables
c = 343;            % Wave speed (m/s)
rho = 1;
h = c * k;          % Grid spacing (m)
L = 3;              % Length

N = floor(L/h);             % Number of points (-)
h = L/N;                    % Recalculate gridspacing from number of points

lambda = c * k / h

% a1 = 1; % loss

% Set cross-sectional geometry
[S, SHalf, SBar] = setTube (N);

%% Lip variables
f0 = 200 / (2*pi);                   % fundamental freq lips
M = 5.37e-5;                  % mass lips
omega0Init = 2 * pi * f0;  % angular freq

%% viscothermal effects
T = 26.85;
[c, rho, eta, nu, gamma] = calcThermoDynConstants(T);
  
sigInit = 5;
H0 = 2.9e-4;

y = 0;
yPrev = y;

w = 1e-2;
Sr = 1.46e-5;

%Initialise states
pNext = zeros(N, 1);
p = zeros(N, 1);
vNext = zeros(N-1, 1);
v = zeros(N-1, 1);

amp = 3000;
% if ~impulse  
%     % input signal
%     t = (0:lengthSound - 1) / fs;
%     freq = 446/4;
%     in = cos(2 * pi * freq * t) - 0.5;
%     in = (in + abs(in)) / 2; % subplus
%     in = in - sum(in) / length(in);
%     in = in * amp;
%     rampLength = 1000; 
%     env = [linspace(0, 1, rampLength), ones(1, lengthSound - rampLength)];
%     in = in .* env;
%     in = in - sum(in) / length(in);
% else
    in = zeros(lengthSound, 1);
%     p(floor(N / 3) - 5 : floor(N / 3) + 5) = 100*hann(11);
% end

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
kinEnergy = zeros (lengthSound, 1);
potEnergy = zeros (lengthSound, 1);
totEnergy = zeros (lengthSound, 1);
qHReed = zeros (lengthSound, 1);
pHReed = zeros (lengthSound, 1);
qHReed = zeros (lengthSound, 1);
pHReed = zeros (lengthSound, 1);

for n = 1:lengthSound
    % if lipstate + equilibrium is below 0, change heaviside variable to 1
    if (y + H0) > 0
        theta = 0;
    else
        theta = 1;
    end
    theta = 0;
    omega0Init = 2 * pi * f0 * (1 + 0.1 * n / fs);
    omega0 = omega0Init * sqrt(1 + 3 * theta);

    sig = sigInit * (1 + 4 * theta);
%     omega0Save(n) = omega0;
    %% Calculate velocities before lip model
    vNext(vRange) = v(vRange) - lambda / (rho * c) * (p(vRange+1) - p(vRange));

    %% Lip model
    
%     if y > 0
%         theta = 0;
%     else
%         theta = 1;
%     end
%     
%     ramp = 1000;
%     if n < ramp
%         Pm = amp * n / ramp;
%     else
        Pm = amp;
%     end
%     
    a1 = 2 / k + omega0^2 * k + sig;
    a2 = Sr / M;
    a3 = 2/k * 1/k * (y - yPrev) - omega0^2 * yPrev;
    b1 = SHalf(1) * vNext(1) + h * SBar(1) / (rho * c^2 * k) * (Pm  - p(1));
    b2 = h * SBar(1) / (rho * c^2 * k);
    c1 = w * subplus(y) * sqrt(2 / rho);
    c2 = b2 + a2 * Sr / a1;
    c3 = b1 - a3 * Sr / a1;
    
    deltaP = sign(c3) * ((-c1 + sqrt(c1^2 + 4 * c2 * abs(c3)))/ (2 * c2))^2;
    
    alpha = 4 / (2 + omega0^2 * k^2 + sig * k);
    beta = (sig * k - 2 - omega0^2 * k^2) / (2 + omega0^2 * k^2 + sig * k);
    epsilon = 2 * Sr * k^2 / (M * (2 + omega0^2 * k^2 + sig * k));
% epsilon = 0;
    yNext(n) = alpha * y + beta * yPrev + epsilon * deltaP;
%     if n < 100
%         yNext(n) = 0;
%     end
    Ub = w * subplus(y + H0) * sign(deltaP) * sqrt(2 * abs(deltaP)/rho);
    Ur = Sr * 1/(2*k) * (yNext(n) - yPrev);
    UbSave(n) = Ub;
    UrSave(n) = Ur;

    %% Schemes

    % calculate schemes
    pNext(pRange) = p(pRange) - rho * c * lambda ./ SBar(pRange) .* (SHalf(pRange) .* vNext(pRange) - SHalf(pRange-1) .* vNext(pRange-1));
    pNext(1) = p(1) - rho * c * lambda ./ SBar(1) .* (-2 * (Ub + Ur) + 2 * SHalf(1) * vNext(1));

%     pNext(N) = p(N) - rho * c * lambda ./ SBar(N) .* (-2 * SHalf(end) .* vNext(end));
    
    % set output from output position
    out(n) = p(outputPos);
    
    %% Energies
    kinEnergy(n) = 1/(2 * rho * c^2) * h * sum(SBar .* scaling .* p.^2);
    potEnergy(n) = rho / 2 * h * sum(SHalf .* vNext .* v);
    hTube(n) = kinEnergy(n) + potEnergy(n);
    hReed(n) = M / 2 * ((1/k * (y - yPrev))^2 + omega0^2 * (y^2 + yPrev^2) / 2);
    qReed(n) = M * sig * (1/(2*k) * (yNext(n) - yPrev))^2 + w * subplus(y + H0) * sqrt(2 / rho) * abs(deltaP)^(3/2);
    idx = n - (1 * (n~=1));
    qHReed(n) = k * qReed(n) + qHReed(idx);
    pReed(n) = -(Ub + Ur) * Pm;
    pHReed(n) = k * pReed(n) + pHReed(idx);
    
    totEnergy(n) = hTube(n) + hReed(n) + qHReed(idx) + pHReed(idx);
    if n>2
        scaledTotEnergy(n) = (totEnergy(n) - hTube(2) - hReed(2)) / (hTube(2) + hReed(2));
    end
    % draw things
    if drawThings && mod (n, drawSpeed) == 0
        subplot(4,1,1)
        cla
        hold on;
%         plotPressurePotential (p / 10000, sqrt(S));
        plot(p)
%         plot(sqrt(S), 'k');
%         plot(-sqrt(S), 'k');
        xlim([1 N]);
        scatter(1, (y + H0) * 10000)
%         ylim([-max(sqrt(S)) max(sqrt(S))] * 1.1);
        title("Pressure");
        
        subplot(4,1,2)
        
        cla;
        plot(vNext);
        hold on;
        plot(sqrt(S) * amp, 'k');
        plot(-sqrt(S) * amp, 'k');
        xlim([1 N]);
        title("Particle Velocity")

        subplot(4,1,3)
        plot(yNext(1:n));
        
        subplot(4,1,4)
%         if n>2
%             plot(scaledTotEnergy(1:n));
%         end
%         title("Normalised total energy (should be 0 within machine precision)")
        hold off;
        plot(hTube(2:n) - hTube(2));
        hold on;
        plot(hReed(2:n) - hReed(2));
        plot(scaledTotEnergy(2:n)*1e-5);

        drawnow;
        
    end
   
    % update states
    v = vNext;
    p = pNext;
    
    yPrev = y;
    y = yNext(n);
end   
plot(yNext)

function [S, SHalf, SBar] = setTube(N)
    mp = linspace(0.0005, 0.0005, floor(N/40));       % mouthpiece
    m2t = linspace(mp(end), 0.0001, floor(N/40));     % mouthpiece to tube
    alpha = 0.25;
    b = m2t(end) * exp(alpha * (0:25));               % bell
    pointsLeft = N - length([mp, m2t, b]);
    tube = linspace(m2t(end), m2t(end), pointsLeft);        % tube

    S = [mp, m2t, tube, b]';                        % True geometry
%     S = 5000 * ones(N,1);
    % Calculate approximations to the geometry
    SHalf = (S(1:N-1) + S(2:N)) * 0.5;           	% mu_{x+}
    SBar = (SHalf(1:end-1) + SHalf(2:end)) * 0.5;
    SBar = [S(1); SBar; S(end)];                    % mu_{x-}S_{l+1/2}
    
end
