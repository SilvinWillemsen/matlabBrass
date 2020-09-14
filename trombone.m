%{
    Full trombone: lip, dynamic tube with proper geometry, radiation 
%}

clear all;
close all;

% drawing variables
drawThings = true;
drawSpeed = 10000;
drawStart = 0;
drawSpeedInit = drawSpeed;
centered = true;

changeL = false;
changeF0 = false;
radiation = true;

connectedToLip = true;

fs = 44100;             % Sample rate (Hz)
k = 1/fs;               % Time step (s)
lengthSound = fs;   % Duration (s)

%% viscothermal effects
T = 26.85;
[c, rho, eta, nu, gamma] = calcThermoDynConstants (T);

%% Tube variables
h = c * k;              % Grid spacing (m)
LnonExtended = 2.658;
%%% Lines added for perfect N (and so perfect energy)
% Ninit = floor(LnonExtended / h);
% LnonExtended = Ninit * h;
%%%
Ninit = LnonExtended / h;
NnonExtended = floor(Ninit);
L = Ninit * h;          % Length
LInit = L;

N = floor(L/h);         % Number of points (-)
alf = Ninit - N;
N = floor(L/h);         % Number of points (-)

lambda = c * k / h      % courant number

%% Lip Collision
Kcol = 0;
alfCol = 5; 

%% Set cross-sectional geometry
[S, SHalf, SBar] = setTube (N+1, NnonExtended, 0);

% Quick note: N is the number of spaces between the points so the number of points is N+1

%% Melody
melody = 0:-1:-12;
melody = [melody, fliplr(melody)];
multiplier = 2.^(-melody./ 12);
range = 1:ceil(lengthSound / length(melody));
multiplierRange = zeros(lengthSound, 1);

% freqs = [300, 283.1623, 267.27, 252.2689];
% lengths = [1.1254, 1.1852, 1.2590, 1.3282];
lengths = LnonExtended * multiplier;
for i = 1:length(melody)
%     freqsRange(range + length(range) * (i-1)) = freqs(i);
    lengthRange(range + length(range) * (i-1)) = lengths(i);
end

%% Lip variables
f0Init = 300;                  % fundamental freq lips
M = 5.37e-5;                % mass lips
omega0 = 2 * pi * f0Init;   % angular freq

sig = 5;                % damping
H0 = 2.9e-4;                % equilibrium
y = 0;                      % initial lip state

if connectedToLip
    w = 1e-2;                   % lip width
    Sr = 1.46e-5;               % lip area
    yPrev = H0;                  % previous lip state
else
    w = 0;
    Sr = 0;         
    yPrev = 0;
end

amp = 3000;                 % input pressure (Pa)

%% Initialise states
upNext = zeros(ceil(N/2) + 1, 1);
up = zeros(ceil(N/2) + 1, 1);
uvNext = zeros(ceil(N/2), 1); 
uv = zeros(ceil(N/2), 1);

if ~connectedToLip
    inputRange = floor(length(up) / 4 - 5):floor(length(up)/4) + 5;
    up(floor(inputRange)) = up(inputRange) + hann(11);
%     up(1:end-1) = rand(length(up)-1, 1);
end

wpNext = zeros(floor(N/2) + 1, 1);
wp = zeros(floor(N/2) + 1, 1);
wvNext = zeros(floor(N/2), 1); % the total v vector is one smaller than the total p vector
wv = zeros(floor(N/2), 1);

% Initialise output
out = zeros (lengthSound, 1);
outputPos = floor(4/5 * N);

% Set ranges
upRange = 2:length(up)-1;         % range without boundaries
wpRange = 2:length(wp)-1;

%% Initialise energies
kinEnergyU = zeros (lengthSound, 1);
potEnergyU = zeros (lengthSound, 1);
hTubeU = zeros (lengthSound, 1);

kinEnergyW = zeros (lengthSound, 1);
potEnergyW = zeros (lengthSound, 1);
hTubeW = zeros (lengthSound, 1);

hReed = zeros (lengthSound, 1);
hColl = zeros (lengthSound, 1);
hRad = zeros (lengthSound, 1);

qReed = zeros (lengthSound, 1);
qHReed = zeros (lengthSound, 1);
pReed = zeros (lengthSound, 1);
pHReed = zeros (lengthSound, 1);
qRad = zeros (lengthSound, 1);
qHRad = zeros (lengthSound, 1);

dampEnergy = zeros (lengthSound, 1);
totH = zeros (lengthSound, 1);
scaledTotEnergy = zeros (lengthSound, 1);

potScalingU = ones(length(up),1);
potScalingU(1) = 1 / 2;
potScalingU(end) = 1 / 2;

potScalingW = ones(length(wp),1);
potScalingW(1) = 1 / 2;
potScalingW(end) = 1 / 2;

ip = [alf * (alf - 1) * (alf - 2) / -6, ...
                (alf - 1) * (alf + 1) * (alf - 2) / 2, ...
                alf * (alf + 1) * (alf - 2) / -2, ...
                alf * (alf + 1) * (alf - 1) / 6];
         
SBarI = SBar(length(uv)-1:length(uv)+2) .* ip';
uvMph = 0;
wvmh = 0;

psiPrev = 0;
etaC = 0;

%% Radiation impedance
R1 = rho * c;
rL = sqrt(SBar(end)) / (2 * pi);
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

filterCoeff = 0.9999;
Pmprev = 0;
flag = false;
% multiplierVec = reshape(repmat(multiplier, lengthSound / 4, 1), lengthSound, 1);
for n = 1:lengthSound
    [S, SHalf, SBar] = setTube (N+1, NnonExtended, n);

    % change wave speed
    if changeL
        LPrev = L;
        filterCoeff = 0.9999;
%         L = (1-filterCoeff) * LInit * 1/multiplierVec(n) + filterCoeff * LPrev;%* sin(0.5 * pi * n/fs));
%         L = LInit * (1 + 0.5 * n / fs);
        L = LInit * (2^sin(2 * pi * n/lengthSound));
%         L = (1-filterCoeff) * lengthRange(n) + filterCoeff * LPrev;
    else
        L = L;
    end
    
    if changeF0
%         f0 = f0Init * multiplierVec(n);
%         f0 = freqsRange(n);
        f0 = f0Init * (1 + 0.05 * sin(2 * pi * 4 * n / fs));
        omega0 = 2 * pi * f0;
    else
        f0 = f0Init;
        omega0 = 2 * pi * f0;

    end
    LSave(n) = L;
    f0Save(n) = f0;
    % save previous state for comparison later
    NPrev = N;

    % recalculate gridspacing, points lambda^2 and alpha from new wave speed
    h = c*k;
    Ninit = L/h;
    N = floor(L/h);
    Nsave(n) = N;
    hSave(n) = h;
    
    lambda = c * k / h;

    alf = (Ninit - N);
%     if (alf < 0.1 && ~flag)
%         drawSpeed = 1;
%         flag = true;
%     elseif (alf > 0.1 && alf < 0.9 && flag)
%         drawSpeed = drawSpeedInit;
%         flag = false;
%     end
    % calculate interpolator
%     if interpol == "cubic"
        ip = [alf * (alf - 1) * (alf - 2) / -6, ...
                (alf - 1) * (alf + 1) * (alf - 2) / 2, ...
                alf * (alf + 1) * (alf - 2) / -2, ...
                alf * (alf + 1) * (alf - 1) / 6];
%     else
%         ip = [0, (1-alf), alf, 0];
%     end

    if abs(N - NPrev) > 1
        disp('too fast')
    end
    
    % add point if N^n > N^{n-1}
    if N > NPrev
        customIp = [-alf * (alf + 1) / ((alf + 2) * (alf + 3)), ...
                    2 * alf * (alf + 1) / ((alf + 1) * (alf + 2)), ...
                    2 * (alf + 1) / ((alf + 2) * (alf + 1)), ...
                    -2 * alf / ((alf + 3) * (alf + 2))];
        if mod(N,2) == 1
%             uvNext = [uvNext; (ip(4) * uvNext(end-1) + ip(3) * uvNext(end) + ip(2) * wvNext(1) + ip(1) * wvNext(2))];
%             uv = [uv; (ip(4) * uv(end-1) + ip(3) * uv(end) + ip(2) * wv(1) + ip(1) * wv(2))];
            uvNext = [uv; uvNextMph];
            uv = [uv; uvMph];
            uvMph =  customIp * [uv(end-1:end); wv(1:2)];
            upNext = [upNext; customIp * [upNext(end-1:end); wpNext(1:2)]];
            up = [up; customIp * [up(end-1:end); wp(1:2)]];
            
        else 
%             wvNext = [(ip(1) * uvNext(end-1) + ip(2) * uvNext(end) + ip(3) * wvNext(1) + ip(4) * wvNext(2)); wvNext];
%             wv = [(ip(1) * uv(end-1) + ip(2) * uv(end) + ip(3) * wv(1) + ip(4) * wv(2)); wv];
            wvNext = [wvNextmh; wvNext];
            wv = [wvmh; wv];
            wvmh =  fliplr(customIp) * [uv(end-1:end); wv(1:2)];
            wpNext = [(customIp(4) * upNext(end-1) + customIp(3) * upNext(end) + customIp(2) * wpNext(1) + customIp(1) * wpNext(2)); wpNext];
            wp = [(customIp(4) * up(end-1) + customIp(3) * up(end) + customIp(2) * wp(1) + customIp(1) * wp(2)); wp];
               
        end
        [S, SHalf, SBar] = setTube(N+1, NnonExtended,n);
        % insert matrix creation here
        
        potScalingU = ones(length(up),1);
        potScalingW = ones(length(wp),1);
        potScalingU(1) = 0.5;
        potScalingW(end) = 0.5;
        potScalingU(end) = 0.5;
        potScalingW(1) = 0.5;
    end   
    
    % remove point if N^n < N^{n-1}
    if N < NPrev
        if mod(N,2) == 0
            uvMph = uv(end);
            uvNext = uvNext(1:end-1);
            uv = uv(1:end-1);
            upNext = upNext(1:end-1);
            up = up(1:end-1);

        else 
            wvmh = wv(1);
            wvNext = wvNext(2:end);
            wv = wv(2:end);
            wpNext = wpNext(2:end);
            wp = wp(2:end);

        end
        if flag
            disp("point removed")
        end
        [S, SHalf, SBar] = setTube(N+1, NnonExtended, n);
        
        potScalingU = ones(length(up),1);
        potScalingW = ones(length(wp),1);
        potScalingU(1) = 0.5;
        potScalingW(end) = 0.5;
        potScalingU(end) = 0.5;
        potScalingW(1) = 0.5;

    end
    upRange = 2:length(up)-1;         % range without boundaries
    wpRange = 2:length(wp)-1;

    A = [1, -ip(4); ...
         -ip(4), 1];
    v = [ip(1) * up(end-2) + ip(2) * up(end-1) + ip(3) * up(end); ...
         ip(3) * wp(1) + ip(2) * wp(2) + ip(1) * wp(3)];
    solut = A \ v;
    
    %% Calculate velocities
    uvNext = uv - lambda / (rho * c) * (up(2:end) - up(1:end-1));
    uvNextMph = uvMph - lambda / (rho * c) * (solut(2) - up(end));
    
    wvNext = wv - lambda / (rho * c) * (wp(2:end) - wp(1:end-1));
    wvNextmh = wvmh - lambda / (rho * c) * (wp(1) - solut(1));
    
    %% Variable input force
%     Pm = filterCoeff * Pmprev + (1 - filterCoeff) * amp;
%     Pmprev = Pm;
%     if n+100 <= lengthSound
% 
%         if lengthRange(n+1) ~= lengthRange(n)
%             if ~connectedToLip
% %                 inputRange = floor(length(up) / 4 - 5):floor(length(up)/4) + 5;
% %                 up(floor(inputRange)) = up(inputRange) + hann(11);
%             else
%                 Pmprev = 0;
%             end
%         end
%     end
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

    %% Calculate pressure
    %% Calculate pressure
    upNext(upRange) = up(upRange) - rho * c * lambda ./ SBar(upRange) .* (SHalf(upRange) .* uvNext(upRange) - SHalf(upRange-1) .* uvNext(upRange-1));
    upNext(1) = up(1) - rho * c * lambda / SBar(1) .* (-2 * (Ub + Ur) + 2 * SHalf(1) * uvNext(1));
    upNext(end) = up(end) - rho * c * lambda ./ SBar(length(up)) .* (SHalf(length(up)) .* uvNextMph - SHalf(length(up) - 1) .* uvNext(end));
    
    wpNext(wpRange) = wp(wpRange) - rho * c * lambda ./ SBar(wpRange + length(up) - 1) .* (SHalf(wpRange + length(up) - 1) .* wvNext(wpRange) - SHalf(wpRange + length(up) - 2) .* wvNext(wpRange-1));
    wpNext(1) = wp(1) - rho * c * lambda ./ SBar(length(up)) .* (SHalf(length(up)) .* wvNext(1) - SHalf(length(up) - 1) .* wvNextmh);
    if radiation
        wpNext(end) = ((1 - rho * c * lambda * z3) * wp(end) - 2 * rho * c * lambda * (v1 + z4 * p1 - (SHalf(end) .* wvNext(end))/SBar(end))) / (1 + rho * c * lambda * z3);
    end
    v1Next = v1 + k / (2 * Lr) * (wpNext(end) + wp(end));
    p1Next = z1 / 2 * (wpNext(end) + wp(end)) + z2 * p1;

    %% Set output from output position
    out(n) = wp(end-3);
    
    %% Energies
    kinEnergyU(n) = rho / 2 * h * sum(SHalf(1:length(uv)) .* uvNext .* uv);
    potEnergyU(n) = 1/(2 * rho * c^2) * h * sum (SBar(1:length(up)) .* potScalingU .* up.^2);
    hTubeU(n) = kinEnergyU(n) + potEnergyU(n);
    
    kinEnergyW(n) = rho / 2 * h * sum(SHalf(length(uv)+1:end) .* wvNext .* wv);
    potEnergyW(n) = 1/(2 * rho * c^2) * h * sum (SBar(length(up):end) .* potScalingW .* wp.^2);
    hTubeW(n) = kinEnergyW(n) + potEnergyW(n);

    hReed(n) = M / 2 * ((1/k * (y - yPrev))^2 + omega0^2 * (y^2 + yPrev^2) / 2);
    hColl(n) = psiPrev^2 / 2;
    hRad(n) = SBar(end) / 2 * (Lr * v1^2 + Cr * p1^2);

    v3Next = p1Next / R2;
    v3 = p1 / R2;
    pBar = 0.5 * (wpNext(end) + wp(end));
    muTPv2 = (pBar - 0.5 * (p1Next + p1)) / R1;
    
    % summed forms (damping and power input)
    idx = n - (1 * (n~=1));
    qReed(n) = M * sig * (1/(2*k) * (yNext(n) - yPrev))^2 + Ub * deltaP;
    qHReed(n) = k * qReed(n) + qHReed(idx);
    pReed(n) = -(Ub + Ur) * Pm;
    pHReed(n) = k * pReed(n) + pHReed(idx);
    qRad(n) = SBar(end) * (R1 * muTPv2^2 + R2 * (0.5 * (v3Next + v3))^2);
    qHRad(n) = k * qRad(n) + qHRad(idx);

    % total energies
    totH(n) = hTubeU(n) + hTubeW(n) + hReed(n) + hColl(n) + hRad(n);
    dampEnergy(n) = qHReed(idx) + pHReed(idx) + qHRad(idx);
    scaledTotEnergy(n) = (totH(n) - totH(1) + dampEnergy(n)) / 2^floor(log2(totH(1)));

    %% Draw things
    if drawThings && mod (n, drawSpeed) == 0 && n> drawStart
        hLocsLeft = (0:(length(up))-1) * h;
        hLocsRight = flip(L - ((0:(length(wp)-1)) * h));   
%         % Plot the velocity
%         subplot(4,1,1)
%         cla
%         hold on;
% %         plotPressurePotential (p / 10000, sqrt(S));
% %         plot(p / 100000)
% %         plot(sqrt(S), 'k');
% %         plot(-sqrt(S), 'k');
%         plot(hLocs * N / L, p, '-o');
% %         xlim([1 N]);
% %         scatter(1, (y + H0) * 10000)
% %         ylim([-max(sqrt(S)) max(sqrt(S))] * 1.1);
% %         title("Pressure");
%         
% %         % Plot the velocity
% %         subplot(4,1,2)
% %         cla;
%         plot(hLocs(1:end-1) * N / L + 0.5, vNext * 100, 'Marker', '.', 'MarkerSize', 10);
% %         hold on;
% %         plot(sqrt(S) * amp, 'k');
% %         plot(-sqrt(S) * amp, 'k');
% %         xlim([1 N]);
% %         title("Particle Velocity")
%         pause(0.2)
%         % Plot the output
% %         subplot(4,1,3)
% %         plot(out(1:n))
% %         
% %         % Plot scaled energy
% %         subplot(4,1,4)
% %         plot(scaledTotEnergy(2:n))
% % %         plot(totEnergy(10:n) - hTube(1) - hReed(1) - hColl(1) - hRad(1))
% Plot the velocity
        subplot(3,1,1)
        cla
        hold on;
        plot(hLocsLeft / L, up, '-o');
        plot(hLocsRight / L, wp, '-o');
        plot(hLocsLeft / L + 0.5 / N, [uvNext; uvNextMph] * 100, 'Marker', '.', 'MarkerSize', 10, 'Color', 'r');
        plot(hLocsRight / L - 0.5 / N, [wvNextmh; wvNext] * 100, 'Marker', '.', 'MarkerSize', 10,  'Color', 'b');
        plot(hLocsLeft(end) / L + 0.5 / N, uvNextMph * 100, 'Marker', 'o', 'MarkerSize', 10, 'Color', 'r');
        plot(hLocsRight(1) / L - 0.5 / N, wvNextmh * 100, 'Marker', 'o', 'MarkerSize', 10,  'Color', 'b');

        %         plot([hLocsLeft, hLocsRight(2:end)] / L, sqrt(S / pi) * 10, 'k');
%         plot([hLocsLeft, hLocsRight(2:end)] / L, -sqrt(S / pi) * 10, 'k');

%         xlim([0.49 0.51])
        ylim([-10000 10000])
%         Plot scaled energy
        subplot(3,1,2)
        plot(out(1:n))
%         hold off
%         plot(kinEnergyU(1:n) + kinEnergyW(1:n))
%         hold on;
%         plot(potEnergyU(1:n) + potEnergyW(1:n))
%         plot(totEnergy(1:n) / totEnergy(1) - 1)
        subplot(3,1,3)
        plot(scaledTotEnergy(2:n))
%         plot(totH(1:n));
        pause(0.5)
        drawnow;
        
        
    end

    %% Update states
    uv = uvNext;
    up = upNext;
        
    wv = wvNext;
    wp = wpNext;

    uvMph = uvNextMph;
    wvmh = wvNextmh;
    
    p1 = p1Next;
    v1 = v1Next;
    
    yPrev = y;
    y = yNext(n);
    psiPrev = psi; 

end

function [S, SHalf, SBar] = setTube(N, NnonExtended, n)
    lengths = [0.708, 0.177, 0.711, 0.241, 0.254, 0.502];
    radii = [0.0069, 0.0072, 0.0069, 0.0071, 0.0075, 0.0107]; % two radii for tuning slide

    lengthN = round(NnonExtended * lengths ./ sum(lengths));
    inner1 = ones(lengthN(1), 1) * radii(1);
    inner2 = ones(lengthN(3), 1) * radii(3);
    gooseneck = ones(lengthN(4), 1) * radii(4);
    tuning = linspace(radii(5), radii(6), lengthN(5))';
    
    x0 = 0.0174; 
    b = 0.0063;
    flare = 0.7;
    bellL = lengthN(end);
    
    bell = b * ((0.502:-0.502 / (bellL - 1):0) + x0).^(-flare);
%     pointsLeft = N - length([mp, m2t, bell]);
%     tube = linspace(m2t(end), m2t(end), pointsLeft);    % tube
    totLengthN = length(inner1) + length(inner2) + length(gooseneck) + length(tuning) + length(bell);

    slide = ones(N - totLengthN, 1) * radii(2);
    totRadii = [inner1; slide; inner2; gooseneck; tuning; bell'];
    
    % True geometry
    S = totRadii.^2 * pi;
%     S = 1 + 1 * 0.5 * (1 + sin(2 * pi * 100 * n / 44100)) * (ones(length(S), 1));

    % Calculate approximations to the geometry
    SHalf = (S(1:N-1) + S(2:N)) * 0.5;                  % mu_{x+}
    SBar = (SHalf(1:end-1) + SHalf(2:end)) * 0.5;
    SBar = [S(1); SBar; S(end)];                        % mu_{x-}S_{l+1/2}
    
end
