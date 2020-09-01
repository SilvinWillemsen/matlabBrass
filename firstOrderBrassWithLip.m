%{
    First order brass
%}

clear all;
close all;

% drawing variables
drawThings = true;
drawSpeed = 5000;
centered = true;

impulse = true;

fs = 44100;         % Sample rate (Hz)
k = 1/fs;           % Time step (s)
lengthSound = fs * 2; % Duration (s)

%% viscothermal effects
T = 26.85;
[c, rho, eta, nu, gamma] = calcThermoDynConstants(T);

%% Tube variables
h = c * k;          % Grid spacing (m)
L = 1;              % Length

N = floor(L/h);             % Number of points (-)
h = L/N;                    % Recalculate gridspacing from number of points

lambda = c * k / h

% a1 = 1; % loss

% lipcollision
Kcol = 1e20;
alf = 1;

% Set cross-sectional geometry
[S, SHalf, SBar] = setTube (N);

%% Lip variables
f0 = 100;                   % fundamental freq lips
M = 5.37e-5;                  % mass lips
omega0Init = 2 * pi * f0;  % angular freq

sigInit = 5;
H0 = 2.9e-4;

y = 0;
yPrev = -H0;

w = 1e-2;
Sr = 1.46e-5;

% w = 0;
% Sr = 0;

%Initialise states
pNext = zeros(N, 1);
p = zeros(N, 1);
vNext = zeros(N-1, 1);
v = zeros(N-1, 1);

amp = 3000;

in = zeros(lengthSound, 1);
% p(floor(2*N / 3) - 5 : floor(2*N / 3) + 5) = hann(11);

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
uBHReed = zeros (lengthSound, 1);
pHReed = zeros (lengthSound, 1);

psiPrev = 0;
eta = 0;

for n = 1:lengthSound
    %% Calculate velocities before lip model
    vNext(vRange) = v(vRange) - lambda / (rho * c) * (p(vRange+1) - p(vRange));
    
    %% Lip model

    % if lipstate + equilibrium is below 0, change heaviside variable to 1
    if (y + H0) > 0
        theta = 0;
    else
        theta = 1;
    end
    theta = 0;
    omega0Init = 2 * pi * f0;% * (1 + 0.1 * n / fs);
    omega0 = omega0Init * sqrt(1 + 3 * theta);

    sig = sigInit * (1 + 4 * theta);

    ramp = 1000;
    if n < ramp
        Pm = amp * n / ramp;
    else
        Pm = amp;
    end
%       

    barr = -H0;
    eta = barr - y;
    g = 0;
    if alf == 1
        if eta > 0
            g = sqrt(Kcol * (alf+1) / 2);
        end
    else
        g = sqrt(Kcol * (alf+1) / 2) * subplus(eta)^((alf - 1)/2);
    end
    
    gSave(n) = g;
    a1 = 2 / k + omega0^2 * k + sig - k/2 * g^2 / M;
    a2 = Sr / M;
    a3 = 2/k * 1/k * (y - yPrev) - omega0^2 * yPrev - psiPrev * g / M;
    b1 = SHalf(1) * vNext(1) + h * SBar(1) / (rho * c^2 * k) * (Pm  - p(1));
    b2 = h * SBar(1) / (rho * c^2 * k);
    c1 = w * subplus(y + H0) * sqrt(2 / rho);
    c2 = b2 + a2 * Sr / a1;
    c3 = b1 - a3 * Sr / a1;
    
    deltaP = sign(c3) * ((-c1 + sqrt(c1^2 + 4 * c2 * abs(c3)))/ (2 * c2))^2;
    
    divTerm = (2 + g^2 * k^2 / (2 * M) + omega0^2 * k^2 + sig * k);
    alpha = 4 / divTerm;
    beta = (sig * k - 2 - omega0^2 * k^2 + g^2 * k^2 / (2 * M)) / divTerm;
    epsilon = 2 * Sr * k^2 / (M * divTerm);
% epsilon = 0;
    yNext(n) = alpha * y + beta * yPrev + epsilon * deltaP + 2 * g * k^2 / (M * divTerm) * psiPrev;

    psi = psiPrev + 0.5 * g * ((barr - yNext(n)) - (barr - yPrev));

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
    hColl(n) = psiPrev^2 / 2;
    qReed(n) = M * sig * (1/(2*k) * (yNext(n) - yPrev))^2 + Ub * deltaP;

    idx = n - (1 * (n~=1));

    qHReed(n) = k * qReed(n) + qHReed(idx);
    pReed(n) = -(Ub + Ur) * Pm;
    pHReed(n) = k * pReed(n) + pHReed(idx);

    
    totEnergy(n) = hTube(n) + hReed(n) + hColl(n) + qHReed(idx) + pHReed(idx);
    scaledTotEnergy(n) = (totEnergy(n) - hTube(1) - hReed(1) - hColl(1)) / 2^floor(log2(hTube(1) + hReed(1) + hColl(1)));
% end
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
        plot(out(1:n))
        subplot(4,1,4)
%         if n>10
        plot(scaledTotEnergy(10:n))
%         hold off;
%         plot(totEnergy(2:n) - totEnergy(2))
%         hold on;
%         plot(-hColl(2:n))
%         plot(totEnergy(2:n) - totEnergy(2))
%         hold on;
%         plot(hColl)
%     
        
%         hold on;
%         plot(uBHReed(1:n-1));
%         end
        drawnow;
        
    end
   
    % update states
    v = vNext;
    p = pNext;
    
    yPrev = y;
    y = yNext(n);
    psiPrev = psi; 

end   
plot(out)

function [S, SHalf, SBar] = setTube(N)
    mp = linspace(0.0005, 0.0005, floor(N/40));       % mouthpiece
    m2t = linspace(mp(end), 0.0001, floor(N/40));     % mouthpiece to tube
    alpha = 0.25;
    b = m2t(end) * exp(alpha * (0:25));               % bell
    pointsLeft = N - length([mp, m2t, b]);
    tube = linspace(m2t(end), m2t(end), pointsLeft);        % tube

    S = [mp, m2t, tube, b]';                        % True geometry
    S = 0.005* ones(N,1);
    % Calculate approximations to the geometry
    SHalf = (S(1:N-1) + S(2:N)) * 0.5;           	% mu_{x+}
    SBar = (SHalf(1:end-1) + SHalf(2:end)) * 0.5;
    SBar = [S(1); SBar; S(end)];                    % mu_{x-}S_{l+1/2}
    
end
