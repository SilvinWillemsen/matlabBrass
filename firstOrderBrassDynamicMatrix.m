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
Ninit = 80;
L = Ninit * h;          % Length
N = floor(L/h);         % Number of points (-)
alf = Ninit - N;

lambda = c * k / h      % courant number

%% Lip Collision
Kcol = 100;
alfCol = 5; 

%% Set cross-sectional geometry
[S, SHalf, SBar] = setTube (N);

%% Lip variables
f0 = 120;                   % fundamental freq lips
M = 5.37e-5;                % mass lips
omega0 = 2 * pi * f0;       % angular freq

sig = 5;                    % damping
H0 = 2.9e-4;                % equilibrium

y = 0;                      % initial lip state
yPrev = 0;                  % previous lip state

w = 0;                   % lip width
Sr = 0;               % lip area

%% Initialise states
upNext = zeros(ceil(N/2) + 1, 1);
up = zeros(ceil(N/2) + 1, 1);
uvNext = zeros(ceil(N/2), 1); 
uv = zeros(ceil(N/2), 1);

up(floor(length(up) / 2 - 5):floor(length(up)/2) + 5) = hann(11);

wpNext = zeros(floor(N/2) + 1, 1);
wp = zeros(floor(N/2) + 1, 1);
wvNext = zeros(floor(N/2) + 1, 1); 
wv = zeros(floor(N/2) + 1, 1);

amp = 30000;                 % input pressure (Pa)

in = zeros(lengthSound, 1);

% Initialise output
out = zeros (lengthSound, 1);

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

%% Create matrices
euv = ones(length(up), 1);
Dxuv = lambda / (rho * c) * spdiags([-euv euv], 0:1, length(uv),length(up));
Dxup = rho * c * lambda * spdiags([-[SHalf(1:length(up)-1)./SBar(2:length(up))] ...
    [2 * SHalf(1) / SBar(1); SHalf(2:length(up)-1)./SBar(2:length(up)-1)]], -1:0, length(up),length(uv));

ewv = ones(length(wp), 1);
lastPTerm = (2 * SHalf(length(up) - 1) / SBar(N)) / (1 + rho * c * lambda * z3);
Dxwv = lambda / (rho * c) * spdiags([-ewv ewv], 0:1, length(wv),length(wp));
Dxwp = rho * c * lambda * spdiags([-[SHalf(length(up)-2:end-1)./SBar(length(up)-1:end-1); 0] ...
    [SHalf(length(up)-1:end)./SBar(length(up)-1:end-1); lastPTerm]], -1:0, length(wp),length(wv));

%% Initialise ranges
pRange = 2:N-1;
vRange = 1:N-1;
for n = 1:lengthSound
    ip = [alf * (alf - 1) * (alf - 2) / -6, ...
                (alf - 1) * (alf + 1) * (alf - 2) / 2, ...
                alf * (alf + 1) * (alf - 2) / -2, ...
                alf * (alf + 1) * (alf - 1) / 6];
            
    vInterpolatedPoints = [1, -ip(4); -ip(4), 1] \ [ip(3) * wv(1) + ip(2) * wv(2) + ip(1) * wv(3);
                            ip(1) * uv(end-2) + ip(2) * uv(end-1) + ip(3) * uv(end)];

    pInterpolatedPoints = [1, -ip(4); -ip(4), 1] \ [ip(3) * wp(1) + ip(2) * wp(2) + ip(1) * wp(3);
                            ip(1) * up(end-2) + ip(2) * up(end-1) + ip(3) * up(end)];

    %% Calculate velocities before lip model
    uvNext = uv - Dxuv * up;
    uvNext(end) = uvNext(end) - lambda / (rho * c) * wp(2);
    
    wvNext = wv - Dxwv * wp;
    wvNext(1) = wvNext(1) - lambda / (rho * c) * up(end-1);

    %% Variable input pressure
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
    b1 = SHalf(1) * uvNext(1) + h * SBar(1) / (rho * c^2 * k) * (Pm  - up(1));
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

    %% Calculate pressure (matrix form)
    upNext = up - Dxup * uvNext;
    upNext(1) = upNext(1) + 2 * rho * c * lambda / SBar(1) * (Ub + Ur);    
    upNext(end) = upNext(end) - rho * c * lambda * SHalf(length(up)) / SBar(length(up)) * wvNext(2);

    wpNext = wp - Dxwp * wvNext;
    wpNext(end) = wpNext(end) - wp(end) + ((1 - rho * c * lambda * z3) * wp(end)...
        - (2 * rho * c * lambda * (v1 + z4 * p1))) / (1 + rho * c * lambda * z3);
    wpNext(1) = wpNext(1) - rho * c * lambda * SHalf(length(up) - 1) / SBar(length(up) - 1) * uvNext(end-1);

    
    % non-matrix form
%     pNext(pRange) = p(pRange) - rho * c * lambda ./ SBar(pRange) .* (SHalf(pRange) .* vNext(pRange) - SHalf(pRange-1) .* vNext(pRange-1));
%     pNext(1) = p(1) - rho * c * lambda / SBar(1) .* (-2 * (Ub + Ur) + 2 * SHalf(1) * vNext(1));
%     pNext(N) = ((1 - rho * c * lambda * z3) * p(N) - 2 * rho * c * lambda * (v1 + z4 * p1 - (SHalf(end) .* vNext(end)) / SBar(N))) / (1 + rho * c * lambda * z3);

    v1Next = v1 + k / (2 * Lr) * (wpNext(end) + wp(end));
    p1Next = z1 / 2 * (wpNext(end) + wp(end)) + z2 * p1;

    %% Set output from output position
    out(n) = wp(end-5);
    
%     %% Energies
%     kinEnergy(n) = 1/(2 * rho * c^2) * h * sum(SBar .* scaling .* up.^2);
%     potEnergy(n) = rho / 2 * h * sum(SHalf .* uvNext .* uv);
%     hTube(n) = kinEnergy(n) + potEnergy(n);
%     hReed(n) = M / 2 * ((1/k * (y - yPrev))^2 + omega0^2 * (y^2 + yPrev^2) / 2);
%     hColl(n) = psiPrev^2 / 2;
%     hRad(n) = SBar(N) / 2 * (Lr * v1^2 + Cr * p1^2);
% 
%     v3Next = p1Next / R2;
%     v3 = p1 / R2;
%     pBar = 0.5 * (upNext(N) + up(N));
%     muTPv2 = (pBar - 0.5 * (p1Next + p1)) / R1;
%     
%     % summed forms (damping and power input)
%     idx = n - (1 * (n~=1));
%     qReed(n) = M * sig * (1/(2*k) * (yNext(n) - yPrev))^2 + Ub * deltaP;
%     qHReed(n) = k * qReed(n) + qHReed(idx);
%     pReed(n) = -(Ub + Ur) * Pm;
%     pHReed(n) = k * pReed(n) + pHReed(idx);
%     qRad(n) = SBar(N) * (R1 * muTPv2^2 + R2 * (0.5 * (v3Next + v3))^2);
%     qHRad(n) = k * qRad(n) + qHRad(idx);
% 
%     % total energies
%     totEnergy(n) = hTube(n) + hReed(n) + hColl(n) + hRad(n) + qHReed(idx) + pHReed(idx) + qHRad(idx);
%     totEnergy1 = hTube(1) + hReed(1) + hColl(1) + hRad(1);
% 
%     scaledTotEnergy(n) = (totEnergy(n) - totEnergy1) / 2^floor(log2(totEnergy1));

    %% Draw things
    if drawThings && mod (n, drawSpeed) == 0
        hLocsLeft = (0:(length(up))-1) * h;
        hLocsRight = flip(L - ((0:(length(wp)-1)) * h));
        % Plot the velocity
%         subplot(4,1,1)
        cla
        hold on;
%         plotPressurePotential (up / 10000, sqrt(S));
        plot(hLocsLeft * N / L, up, 'o-')
        plot(hLocsRight * N / L, wp, 'o-')

%         plot(sqrt(S), 'k');
%         plot(-sqrt(S), 'k');
        xlim([1 N]);
%         scatter(1, (y + H0) * 10000)
%         ylim([-max(sqrt(S)) max(sqrt(S))] * 1.1);
%         title("Pressure");
        
        % Plot the velocity
%         subplot(4,1,2)
%         cla;
        plot(hLocsLeft(1:end-1) * N / L + 0.5, uvNext * 100, 'Marker', '.', 'MarkerSize', 10);
        hold on;
        plot(hLocsRight(1:end) * N / L - 0.5, wvNext * 100, 'Marker', '.', 'MarkerSize', 10);
%         pause(0.2);
        if n == 17
            pause;
        end
%         plot(sqrt(S) * amp, 'k');
%         plot(-sqrt(S) * amp, 'k');
%         xlim([1 N]);
%         title("Particle Velocity")
% 
%         % Plot the output
%         subplot(4,1,3)
%         plot(out(1:n))
%         
%         % Plot scaled energy
%         subplot(4,1,4)
% %         plot(scaledTotEnergy(2:n))
% %         plot(totEnergy(10:n) - hTube(1) - hReed(1) - hColl(1) - hRad(1))
        drawnow;
        
    end

    %% Update states
    uv = uvNext;
    up = upNext;
    
    wv = wvNext;
    wp = wpNext;

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
%     S = (ones(length(S), 1));
    % Calculate approximations to the geometry
    SHalf = (S(1:N-1) + S(2:N)) * 0.5;                  % mu_{x+}
    SBar = (SHalf(1:end-1) + SHalf(2:end)) * 0.5;
    SBar = [S(1); SBar; S(end)];                        % mu_{x-}S_{l+1/2}
    
end
