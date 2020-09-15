%{
    First order brass with coupled lip model
%}

clear all;
close all;

% drawing variables
drawThings = true;
drawSpeed = 1;
centered = true;

fs = 44100;             % Sample rate (Hz)
k = 1/fs;               % Time step (s)
lengthSound = fs * 2;   % Duration (s)

%% viscothermal effects
T = 26.85;
[c, rho, eta, nu, gamma] = calcThermoDynConstants (T);

%% Tube variables
h = c * k;              % Grid spacing (m)
L = 0.5;                  % Length

N = floor(L/h);         % Number of points (-)
h = L/N;                % Recalculate gridspacing from number of points

lambda = c * k / h      % courant number

%% Lip Collision
Kcol = 10000;
alfCol = 5; 

%% Set cross-sectional geometry
[S, SHalf, SBar] = setTube (N);

%% Lip variables
f0 = 300;                   % fundamental freq lips
M = 5.37e-5;                % mass lips
omega0 = 2 * pi * f0;   % angular freq

sig = 5;                % damping
H0 = 2.9e-4;                % equilibrium

y = 0;                      % initial lip state
yPrev = H0;                  % previous lip state

w = 1e-2;                   % lip width
Sr = 1.46e-5;               % lip area

%% Initialise states
pNext = zeros(N, 1);        % pressure
p = zeros(N, 1);
% p(N/2-5 : N/2+5) =100 * hann(11);
vNext = zeros(N-1, 1);      % velocity
v = zeros(N-1, 1);

amp = 30000;                 % input pressure (Pa)

in = zeros(lengthSound, 1);
p(floor(2*N / 3) - 5 : floor(2*N / 3) + 5) = 10000 * hann(11);

% Initialise output
out = zeros (lengthSound, 1);
outputPos = floor(4/5 * N);

% Set ranges
pRange = 2:N-1;         % range without boundaries
vRange = 1:N-1;         % range from 1/2 - N-1/2

% Virtual points
SNph = 2 * SBar(N) - SHalf(end);
SOnemh = 2 * SBar(1) - SHalf(1);

%% Initialise energies
kinEnergy = zeros (lengthSound, 1);
potEnergy = zeros (lengthSound, 1);
totEnergy = zeros (lengthSound, 1);
hTube = zeros (lengthSound, 1);
hReed = zeros (lengthSound, 1);
hColl = zeros (lengthSound, 1);
hRad = zeros (lengthSound, 1);

qReed = zeros (lengthSound, 1);
qHReed = zeros (lengthSound, 1);
pReed = zeros (lengthSound, 1);
pHReed = zeros (lengthSound, 1);
qRad = zeros (lengthSound, 1);
qHRad = zeros (lengthSound, 1);


totEnergy = zeros (lengthSound, 1);
scaledTotEnergy = zeros (lengthSound, 1);

scaling = ones(N,1);
if centered
    scaling(1) = 1 / 2;
    scaling(N) = 1 / 2;
end

%
psiPrev = 0;
etaC = 0;

%% Radiation impedance
R1 = rho * c;
rL = sqrt(SBar(N)) / (2 * pi);
Lr = 0.613 * rho * rL;
R2 = 0.505 * rho * c;
Cr = 1.111 * rL / (rho * c^2); 

zDiv = 2 * R1 * R2 * Cr + k * (R1 + R2);
if zDiv == 0
    z1 = 0;
    z2 = 0;
else
    z1 = 2 * R2 * k / zDiv;
    z2 = (2 * R1 * R2 * Cr - k * (R1 + R2)) / zDiv;

end
z3 = k/(2*Lr) + z1 / (2 * R2) + Cr * z1 / k;
z4 = (z2 + 1) / (2 * R2) + (Cr * z2 - Cr) / k;
    
p1 = 0;
v1 = 0;
for n = 1:lengthSound
    %% Calculate velocities before lip model
    vNext = v - lambda / (rho * c) * (p(vRange+1) - p(vRange));
    
    %% Variable input force
    ramp = 1000;
    if n < ramp
        Pm = amp * n / ramp;
    else
        Pm = amp;
    end

    %% Collision
    barr = -H0;
    etaC = barr - y;
    g = 0;
    if alfCol == 1
        if etaC > 0
            g = sqrt(Kcol * (alfCol+1) / 2);
        end
    else
        g = sqrt(Kcol * (alfCol+1) / 2) * subplus(etaC)^((alfCol - 1)/2);
    end
    
    %% Obtain deltaP
    a1 = 2 / k + omega0^2 * k + sig + g^2 * k / (2 * M);
    a2 = Sr / M;
    a3 = 2/k * 1/k * (y - yPrev) - omega0^2 * yPrev + g / M * psiPrev;
    b1 = SHalf(1) * vNext(1) + h * SBar(1) / (rho * c^2 * k) * (Pm  - p(1));
    b2 = h * SBar(1) / (rho * c^2 * k);
    c1 = w * subplus(y + H0) * sqrt(2 / rho);
    c2 = b2 + a2 * Sr / a1;
    c3 = b1 - a3 * Sr / a1;
    
    deltaP = sign(c3) * ((-c1 + sqrt(c1^2 + 4 * c2 * abs(c3)))/ (2 * c2))^2;
    
    %% Update lip scheme
    gammaR = g * k^2 / (2 * M);
    alpha = 2 + omega0^2 * k^2 + sig * k + g * gammaR;
    beta = sig * k - 2 - omega0^2 * k^2 + g * gammaR;
    xi = 2 * Sr * k^2 / M;

    yNext(n) = 4 / alpha * y + beta / alpha * yPrev + xi / alpha * deltaP + 4 * gammaR * psiPrev / alpha;

    %% Update collision potential
    psi = psiPrev - 0.5 * g * (yNext(n) - yPrev);
    
    %% Calculate flow velocities
    Ub = w * subplus(y + H0) * sign(deltaP) * sqrt(2 * abs(deltaP)/rho);
    Ur = Sr * 1/(2*k) * (yNext(n) - yPrev);

    %% Calculate pressure
    pNext(pRange) = p(pRange) - rho * c * lambda ./ SBar(pRange) .* (SHalf(pRange) .* vNext(pRange) - SHalf(pRange-1) .* vNext(pRange-1));
    pNext(1) = p(1) - rho * c * lambda / SBar(1) .* (-2 * (Ub + Ur) + 2 * SHalf(1) * vNext(1));
    pNext(N) = ((1 - rho * c * lambda * z3) * p(N) - 2 * rho * c * lambda * (v1 + z4 * p1 - (SHalf(end) .* vNext(end))/SBar(N))) / (1 + rho * c * lambda * z3);
%     pNext(N) = alphaR * p(N) + betaR * SHalf(end) * vNext(end) + epsilonR * v1 + nuR * p1;

    v1Next = v1 + k / (2 * Lr) * (pNext(N) + p(N));
    p1Next = z1 / 2 * (pNext(N) + p(N)) + z2 * p1;
%     p1Next = taoR * p1 + xiR / 2 * (pNext(N) + p(N)) ;

    %% Set output from output position
    out(n) = p(outputPos);
    
    %% Energies
    kinEnergy(n) = 1/(2 * rho * c^2) * h * sum(SBar .* scaling .* p.^2);
    potEnergy(n) = rho / 2 * h * sum(SHalf .* vNext .* v);
    hTube(n) = kinEnergy(n) + potEnergy(n);
    hReed(n) = M / 2 * ((1/k * (y - yPrev))^2 + omega0^2 * (y^2 + yPrev^2) / 2);
    hColl(n) = psiPrev^2 / 2;
    hRad(n) = SBar(N) / 2 * (Lr * v1^2 + Cr * p1^2);

    v3Next = p1Next / R2;
    v3 = p1 / R2;
    pBar = 0.5 * (pNext(N) + p(N));
    muTPv2 = (pBar - 0.5 * (p1Next + p1)) / R1;
    
    % summed forms (damping and power input)
    idx = n - (1 * (n~=1));
    qReed(n) = M * sig * (1/(2*k) * (yNext(n) - yPrev))^2 + Ub * deltaP;
    qHReed(n) = k * qReed(n) + qHReed(idx);
    pReed(n) = -(Ub + Ur) * Pm;
    pHReed(n) = k * pReed(n) + pHReed(idx);
    qRad(n) = SBar(N) * (R1 * muTPv2^2 + R2 * (0.5 * (v3Next + v3))^2);
    qHRad(n) = k * qRad(n) + qHRad(idx);
%     qRad(n
%     qRad(n) = SBar(N) * (R1 / 2 * (v))

    % total energies
    totEnergy(n) = hTube(n) + hReed(n) + hColl(n) + hRad(n) + qHReed(idx) + pHReed(idx) + qHRad(idx);
    scaledTotEnergy(n) = (totEnergy(n) - hTube(1) - hReed(1) - hColl(1) - hRad(1)) / 2^floor(log2(hTube(1) + hReed(1) + hColl(1) + hRad(1)));

    %% Draw things
    if drawThings && mod (n, drawSpeed) == 0
        hLocs = (0:length(p)-1) * h;
        % Plot the velocity
%         subplot(4,1,1)
        cla
        hold on;
%         plotPressurePotential (p / 10000, sqrt(S));
%         plot(p / 100000)
%         plot(sqrt(S), 'k');
%         plot(-sqrt(S), 'k');
        plot(hLocs * N / L, p, '-o');
%         xlim([1 N]);
%         scatter(1, (y + H0) * 10000)
%         ylim([-max(sqrt(S)) max(sqrt(S))] * 1.1);
%         title("Pressure");
        
%         % Plot the velocity
%         subplot(4,1,2)
%         cla;
        plot(hLocs(1:end-1) * N / L + 0.5, vNext * 100, 'Marker', '.', 'MarkerSize', 10);
%         hold on;
%         plot(sqrt(S) * amp, 'k');
%         plot(-sqrt(S) * amp, 'k');
%         xlim([1 N]);
%         title("Particle Velocity")
pause(0.2)
        % Plot the output
%         subplot(4,1,3)
%         plot(out(1:n))
%         
%         % Plot scaled energy
%         subplot(4,1,4)
%         plot(scaledTotEnergy(2:n))
% %         plot(totEnergy(10:n) - hTube(1) - hReed(1) - hColl(1) - hRad(1))
        drawnow;
        
    end

    %% Update states
    v = vNext;
    p = pNext;
    
    p1 = p1Next;
    v1 = v1Next;
    
    yPrev = y;
    y = yNext(n);
    psiPrev = psi; 

end   

function [S, SHalf, SBar] = setTube(N)
    mp = linspace(0.001, 0.001, floor(N/40));           % mouthpiece
    m2t = linspace(mp(end), 0.0005, floor(N/40));       % mouthpiece to tube
    
    alpha = 0.25;                                       % bell
    b = m2t(end) * exp(alpha * (0:18));
    pointsLeft = N - length([mp, m2t, b]);
    tube = linspace(m2t(end), m2t(end), pointsLeft);    % tube

    S = [mp, m2t, tube, b]';                            % True geometry
%     S = 0.005 * (ones(length(S), 1));
    % Calculate approximations to the geometry
    SHalf = (S(1:N-1) + S(2:N)) * 0.5;                  % mu_{x+}
    SBar = (SHalf(1:end-1) + SHalf(2:end)) * 0.5;
    SBar = [S(1); SBar; S(end)];                        % mu_{x-}S_{l+1/2}
    
end
