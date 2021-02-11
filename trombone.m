%{
    Full trombone: lip, dynamic tube with proper geometry, radiation 
%}

% clear all;
close all;

% drawing variables
drawThings = true;
drawSetting = 1; % change to 0 to go back to previous drawing

drawSpeed = 50;
drawStart = 0;
drawSpeedInit = drawSpeed;

fixedNonInterpolatedL = true;

centered = true;

changeL = ~fixedNonInterpolatedL;
changeF0 = true;
radiation = true;

connectedToLip = true;

fs = 44100;             % Sample rate (Hz)
k = 1/fs;               % Time step (s)
lengthSound = fs * 5;   % Duration (s)

LnonExtended = 2.658;
Lextended = 3.718;

% LnonExtended = 2.685;

% LnonExtended = 0.999959410430839;
%% viscothermal effects
T = 26.85;
[c, rho, eta, nu, gamma] = calcThermoDynConstants (T);

%% Melody
% melody = [0, -1, -3, -5, -7, -8, -10, -12];
melody = [-10, -11, -12, -13, -13];
melody = 12 * log2(300/520);
% melody = [melody, fliplr(melody)];
% multiplier = 2.^(melody./ 12);
multiplier = 1;
range = 1:ceil(lengthSound / length(melody));
multiplierRange = zeros(lengthSound, 1);

freqs = 300;
% freqs = 520 * multiplier;

%% Tube variables
h = c * k;              % Grid spacing (m)
%%% Lines added for perfect N (and so perfect energy)
% Ninit = floor(LnonExtended / h);
% LnonExtended = Ninit * h;
%%%
Ninit = LnonExtended / h;
if fixedNonInterpolatedL
    LnonExtended = floor(Ninit) * h;
    Ninit = LnonExtended / h;
end
% lengths = LnonExtended ./ multiplier;
lengths = LnonExtended ./ multiplier;

NnonExtended = floor(Ninit);

lambda = c * k / h      % courant number
pitchGlide = ones(length(range), 1);
% pitchGlide(1:floor(length(range)*0.3)) = linspace(1.2, 1, length(range) * 0.3);

%% Melody 2
for i = 1:length(melody)
    if i == length(melody)
        pitchGlide = ones(length(range),1);
    end
    freqsRange(range + length(range) * (i-1)) = (2-pitchGlide) .* freqs(i);
    lengthRange(range + length(range) * (i-1)) = pitchGlide .* lengths(i);
end
% L = lengthRange(1);          % Length
L = LnonExtended;
LInit = L;
N = floor(L/h);         % Number of points (-)
alf = Ninit - N;

%% Lip Collision
Kcol = 10000;
alfCol = 3; 

%% Set cross-sectional geometry
[S, SHalf, SBar, addPointsAt] = setTube (N+1, NnonExtended, 0);

% Quick note: N is the number of spaces between the points so the number of points is N+1

%% Lip variables
f0Init = freqs(1);                  % fundamental freq lips
f0 = f0Init;

Mlip = 5.37e-5;                % mass lips
omega0 = 2 * pi * f0Init;   % angular freq

sig = 5;                % damping
H0 = 2.9e-4;                % equilibrium
y = 0;                      % initial lip state

if connectedToLip
%     w = 0.25e-2;                   % lip width
    Sr = 1.46e-5;               % lip area
    w = 0.01;                   % lip width
    yPrev = H0;                  % previous lip state
else
    w = 0;
    Sr = 0;         
    yPrev = 0;
end

if connectedToLip
    amp = 3000;                 % input pressure (Pa)
else
    amp = 100;
end

%% Initialise states
upNext = zeros(ceil(addPointsAt) + 1, 1);
up = zeros(ceil(addPointsAt) + 1, 1);
uvNext = zeros(ceil(addPointsAt), 1); 
uv = zeros(ceil(addPointsAt), 1);

if ~connectedToLip
%     inputRange = floor(length(up) / 4 - 5):floor(length(up)/4) + 5;
    inputRange = 21:31;
    up(floor(inputRange)) = up(floor(inputRange)) + 1000 * (1.0 - cos (2.0 * pi * (0:length(inputRange)-1)' / (length(inputRange) - 1))) * 0.5;
%     up(floor(inputRange)) = up(inputRange) + 500 * hann(11);
%     up(1:end-1) = rand(length(up)-1, 1);
end

wpNext = zeros(floor(N-addPointsAt) + 1, 1);
wp = zeros(floor(N-addPointsAt) + 1, 1);
wvNext = zeros(floor(N-addPointsAt), 1); % the total v vector is one smaller than the total p vector
wv = zeros(floor(N-addPointsAt), 1);

% Initialise output
out = zeros (lengthSound, 1);

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

Pmprev = 0;
flag = false;
statesSave = [];
% multiplierVec = reshape(repmat(multiplier, lengthSound / 4, 1), lengthSound, 1);
for n = 1:lengthSound
    [S, SHalf, SBar] = setTube (N+1, NnonExtended, n);

    filterCoeff = 0.9999;

    if changeL
        LPrev = L;
%         L = (1-filterCoeff) * LInit * 1/multiplierVec(n) + filterCoeff * LPrev;%* sin(0.5 * pi * n/fs));
%         L = LInit * (1 + 0.5 * n / fs);
%         L = LInit * (2^sin(2 * pi * n/lengthSound));
        L = (1-filterCoeff) * lengthRange(n) + filterCoeff * LPrev;
    else
        L = L;
    end
    
    if changeF0
%         f0 = f0Init * multiplierVec(n);
        f0Prev = f0;
        f0 = (1-filterCoeff) * freqsRange(n) + filterCoeff * f0Prev;

%         f0 = freqsRange(n);
%         f0 = f0Init * (1 + 0.05 * sin(2 * pi * 4 * n / fs));
        if n > 3/5 * lengthSound
            omega0 = 2 * pi * f0 * (1 + 0.02 * sin(2 * pi * 3 * n / fs));
        else
            omega0 = 2 * pi * f0;
        end
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
        
        statesSave = [statesSave; [up(end-1), up(end), wp(1), wp(2), uv(end-1), uv(end), wv(1), wv(2), uvMph, wvmh] ];
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
        statesSave = [statesSave; [up(end-1), up(end), wp(1), wp(2), uv(end-1), uv(end), wv(1), wv(2), uvMph, wvmh] ];

    end
    upRange = 2:length(up)-1;         % range without boundaries
    wpRange = 2:length(wp)-1;

    A = [1, -ip(4); ...
         -ip(4), 1];
    v = [ip(1) * up(end-2) + ip(2) * up(end-1) + ip(3) * up(end); ...
         ip(3) * wp(1) + ip(2) * wp(2) + ip(1) * wp(3)];
    solut = A \ v;
    
    quadIp = [-(alf - 1) / (alf + 1), 1, (alf - 1) / (alf + 1)];
    solut = [up(end-1) * quadIp(1) + up(end) * quadIp(2) + wp(1) * quadIp(3); ...
            up(end) * quadIp(3) + wp(1) * quadIp(2) + wp(2) * quadIp(1)];
    
    %% Calculate velocities
    uvNext = uv - lambda / (rho * c) * (up(2:end) - up(1:end-1));
    uvNextMph = uvMph - lambda / (rho * c) * (solut(2) - up(end));
    
    wvNext = wv - lambda / (rho * c) * (wp(2:end) - wp(1:end-1));
    wvNextmh = wvmh - lambda / (rho * c) * (wp(1) - solut(1));
    
%     %% Variable input force
%     filterCoeffPm = 0.9995;
%     Pm = filterCoeffPm * Pmprev + (1 - filterCoeffPm) * amp;
%     Pmprev = Pm;
    Pm = amp * 6;
%     if lengthRange(n+1) ~= lengthRange(n)
    if mod(n, lengthSound / length(melody)) == 0
        if ~connectedToLip
            inputRange = floor(length(up) / 4 - 5):floor(length(up)/4) + 5;
            up(floor(inputRange)) = up(inputRange) + hann(11);
        else
            if n < lengthSound * 4/5 
                Pmprev = 0;
            end
        end
    end
%     ramp = 1000;
%     if n < ramp
%         Pm = amp * n / ramp;
%     else
        Pm = amp;
%     end

    %% Collision
    barr = -H0;
    etaC = barr - y;  
    etaCPrev = barr - yPrev;

    g = 0;
    calcLipDisp; % calculate yNext without collision
    etaCNext = barr - yNext(n);

    if psiPrev < 0
        kappaLip = -1;
    else
        kappaLip = 1;
    end
    
    if etaC >= 0
        g = kappaLip * sqrt(Kcol * (alfCol+1) / 2) * subplus (etaC)^((alfCol - 1.0) / 2.0);
    else
        if(etaCNext - etaCPrev ~= 0)
            g = -2 * psiPrev / (etaCNext - etaCPrev);
        else 
            disp("DIVISION BY 0");
        end
    end
    
    calcLipDisp; % calculate yNext with collision

    %% Update collision potential
    psi = psiPrev - 0.5 * g * (yNext(n) - yPrev);
    
    %% Calculate flow velocities
    Ub = w * subplus(y + H0) * sign(deltaP) * sqrt(2 * abs(deltaP)/rho);
    Ur = Sr * 1/(2*k) * (yNext(n) - yPrev);

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

    if v1Next ~= 0
        disp("wait")
    end
    %% Set output from output position
    out(n) = wp(end);
    
    %% Energies
    kinEnergyU(n) = rho / 2 * h * sum(SHalf(1:length(uv)) .* uvNext .* uv);
    potEnergyU(n) = 1/(2 * rho * c^2) * h * sum (SBar(1:length(up)) .* potScalingU .* up.^2);
    hTubeU(n) = kinEnergyU(n) + potEnergyU(n);
    
    kinEnergyW(n) = rho / 2 * h * sum(SHalf(length(uv)+1:end) .* wvNext .* wv);
    potEnergyW(n) = 1/(2 * rho * c^2) * h * sum (SBar(length(up):end) .* potScalingW .* wp.^2);
    hTubeW(n) = kinEnergyW(n) + potEnergyW(n);

    hReed(n) = Mlip / 2 * ((1/k * (y - yPrev))^2 + omega0^2 * (y^2 + yPrev^2) / 2);
    hColl(n) = psiPrev^2 / 2;
    hRad(n) = SBar(end) / 2 * (Lr * v1^2 + Cr * p1^2);

    v3Next = p1Next / R2;
    v3 = p1 / R2;
    pBar = 0.5 * (wpNext(end) + wp(end));
    muTPv2 = (pBar - 0.5 * (p1Next + p1)) / R1;
    
    % summed forms (damping and power input)
    idx = n - (1 * (n~=1));
    qReed(n) = Mlip * sig * (1/(2*k) * (yNext(n) - yPrev))^2 + Ub * deltaP;
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
        if drawSetting == 0
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
            plot(hLocsLeft / L + 0.5 / N, [uvNext; uvNextMph] / amp, 'Marker', '.', 'MarkerSize', 10, 'Color', 'r');
            plot(hLocsRight / L - 0.5 / N, [wvNextmh; wvNext] / amp, 'Marker', '.', 'MarkerSize', 10,  'Color', 'b');
            plot(hLocsLeft(end) / L + 0.5 / N, uvNextMph / amp, 'Marker', 'o', 'MarkerSize', 10, 'Color', 'r');
            plot(hLocsRight(1) / L - 0.5 / N, wvNextmh / amp, 'Marker', 'o', 'MarkerSize', 10,  'Color', 'b');
%             test1 = find(S(1:length(S) - 66) ~= S(2:length(S)- 65));
%             plot([test1(1) / N; test1(1) / N], 10 * [-amp, amp])
%             plot([test1(2) / N; test1(2) / N], 10 * [-amp, amp])
%             plot([hLocsLeft(1:end-1), hLocsRight] / L, sqrt(S / pi) * 100 * amp, 'k');
%             plot([hLocsLeft(1:end-1), hLocsRight] / L, -sqrt(S / pi) * 100 * amp, 'k');
            %         plot([hLocsLeft, hLocsRight(2:end)] / L, sqrt(S / pi) * 10, 'k');
    %         plot([hLocsLeft, hLocsRight(2:end)] / L, -sqrt(S / pi) * 10, 'k');

    %         xlim([0.49 0.51])
%             ylim([-amp * 10 amp * 10])
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
        else
            subplot(2,1,1)
            hold off;
            plot(1:length(up), up);
%             plot(1:length(uv), uv);
            hold on
            plot(1:M(n)+1, pState(n, 1:M(n)+1), 'Marker', '.', 'MarkerSize', 10);
            plot((1:length(wp)) + length(up) - 1 + alf, wp);
            plot((M(n)+1:(M(n)+Mw(n) + 1)) + alfSave(n), pState(n, (maxM+2):(maxM+Mw(n)+2)), 'Marker', 'o', 'MarkerSize', 2);
%             xlim([M(n) - 5, (M(n)+5)])
%             plot (sqrt(S / pi) * 100 * amp, 'k')
%             hold on
%             plot (sqrt(Ssave(n, :) / pi) * 100 * amp, '--k')
%             hold off
%             plot((1:length(wp)) + length(up) - 1, pState(n, (1:length(wp)) + length(up))');
%             plot((1:length(wv)) + length(uv) - 1 + alf, wv);
% 
            subplot(2,1,2)
            plot([1:length(up), (1:length(wp)) + length(up) - 1 + alf], [up', wp'] - [pState(n, 1:M(n)+1), pState(n, (maxM+2):(maxM+Mw(n)+2)) ])
%             hold off;
%             plot(1:length(up), up, '-o');
%             hold on;   
%             plot((1:length(wp))+length(up)-1, wp, '-o');  
            pause(0.5);
            drawnow;
        end
        
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

function [S, SHalf, SBar, addPointsAt] = setTube(N, NnonExtended, n)
    lengths = [0.708, 0.177, 0.711, 0.306, 0.254, 0.502];
    radii = [0.0069, 0.0072, 0.0069, 0.0071, 0.0075, 0.0107]; % two radii for tuning slide

    lengthN = round(NnonExtended * lengths ./ sum(lengths));
    addPointsAt = round(lengthN(1) + lengthN(2) * 0.5) + (N-NnonExtended) * 0.5; % indicate split of two connected schemes (including offset if N differs from NnonExtended
    
    inner1 = ones(lengthN(1), 1) * radii(1);
    inner2 = ones(lengthN(3), 1) * radii(3);
    gooseneck = ones(lengthN(4), 1) * radii(4);
    tuning = linspace(radii(5), radii(6), lengthN(5))';
    
    x0 = 0.0174; 
    b = 0.0063;
    flare = 0.7;
    bellL = lengthN(end);
    
    bell = b * ((lengths(6):-lengths(6) / (bellL - 1):0) + x0).^(-flare);

    
%     pointsLeft = N - length([mp, m2t, bell]);
%     tube = linspace(m2t(end), m2t(end), pointsLeft);    % tube
    totLengthN = length(inner1) + length(inner2) + length(gooseneck) + length(tuning) + length(bell);
%     lengthN(2) = lengthN(2) + (N - NnonExtended + 1);
    slide = ones(N - totLengthN, 1) * radii(2);
    totRadii = [inner1; slide; inner2; gooseneck; tuning; bell'];
    
    % True geometry
    S = totRadii.^2 * pi;
%     S = 1 + 1 * 0.5 * (1 + sin(2 * pi * 100 * n / 44100)) * (ones(length(S), 1));
%     S = ones(size(S));
    % Calculate approximations to the geometry
    SHalf = (S(1:N-1) + S(2:N)) * 0.5;                  % mu_{x+}
    SBar = (SHalf(1:end-1) + SHalf(2:end)) * 0.5;
    SBar = [S(1); SBar; S(end)];                        % mu_{x-}S_{l+1/2}
    
end
