%{
    Full trombone: lip, dynamic tube with proper geometry, radiation
%}

clear all;
close all;

fs = 44100;             % Sample rate (Hz)
k = 1/fs;               % Time step (s)
lengthSound = fs;

% drawing variables
drawThings = false;
drawSpeed = 10;
drawStart = 0;
drawSpeedInit = drawSpeed;
zoomPlot = false;

shouldDispCorr = false;
correctV = false;
correctVirtual = false;

radiation = true;

onlyAnalysis = true;

LnonExtended = 2.593;
Lextended = 3.653;

Ndiffmax = 5;

%% viscothermal effects
T = 26.85;
[c, rho, eta, nu, gamma] = calcThermoDynConstants (T);

%% Tube variables
h = c * k;              % Grid spacing (m)
NnonExtended = LnonExtended / h;
Nextended = Lextended / h;

Nstart = NnonExtended;
Nend = Nextended;
% Nstart = 150;
% Nend = 21;
if onlyAnalysis
    numLoops = Ndiffmax * abs(Nend - Nstart);
    fSave = zeros(numLoops, ceil(max(Nstart, Nend)));
    sigmaSave = zeros(numLoops, ceil(max(Nstart, Nend)));
else
    numLoops = lengthSound;
end

Ninit = Nstart;

lambdaFact = 1; % actually should be 0.999!
lambda = lambdaFact * c * k / h      % courant number

LInit = Nstart*h;
Lend = Nend * h;

L = LInit;
Ninit = L/h;
N = floor(Ninit);         % Number of points (-)
alf = Ninit - N;

%% Set cross-sectional geometry
setToOnes = false;
[S, SHalf, SBar, addPointsAt] = setTube (N+1, NnonExtended, 0, setToOnes);

% Quick note: N is the number of spaces between the points so the number of points is N+1

%% Initialise states
Mp = ceil(addPointsAt);
Mq = floor(N-addPointsAt);

totSize = Mp+1+Mq+1;

upNext = zeros(Mp + 1, 1);
up = zeros(Mp + 1, 1);
upPrev = zeros(Mp + 1, 1);

wpNext = zeros(Mq + 1, 1);
wp = zeros(Mq + 1, 1);
wpPrev = zeros(Mq + 1, 1);

uvNext = zeros(Mp+1, 1);
uv = zeros(Mp+1, 1);
uvPrev = zeros(Mp+1, 1);

wvNext = zeros(Mq+1, 1);
wv = zeros(Mq+1, 1);
wvPrev = zeros(Mq+1, 1);

pNext = [upNext; wpNext];
p = [up; wp];
pPrev = [upPrev; wpPrev];

vNext = [uvNext; wvNext];
v = [uv; wv];
vPrev = [uvPrev; wvPrev];

pvNext = [pNext; vNext];
pv = [p; v];

uvMphPrev = 0;
wvmhPrev = 0;

upMp1Prev = 0;
wpm1Prev = 0;

exciteU = true;
if exciteU
    width = floor(0.2 * length(up));
    loc = floor(0.5 * length(up));
    inputRange = (loc-width) : (loc+width);
    up(inputRange) = up(inputRange) + hann(width*2+1);
    upPrev = up;
    
    p(inputRange) = p(inputRange) + hann(width*2+1);
    pPrev = p;
    
    pv(inputRange) = pv(inputRange) + hann(width*2+1);

    
else
    width = floor(0.2 * length(wp));
    loc = floor(0.5 * length(wp));
    inputRange = (loc-width) : (loc+width);
    wp(inputRange) = wp(inputRange) + hann(width*2+1);
    wpPrev = wp;

    p(inputRange+Mp+1) = p(inputRange+Mp+1) + hann(width*2+1);
    pPrev = p;
    
    pv(inputRange+1+Mp) = pv(inputRange+1+Mp) + hann(width*2+1);

end

% Initialise output
out = zeros (lengthSound, 1);

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
p1Prev = 0;
v1Prev = 0;

Pmprev = 0;
flag = false;
statesSave = [];
changeL = true;
alternatePV = false;

for n = 1:ceil(numLoops)
    
    retrieveAndRecalculateParams;
    addAndRemovePoints;
    
    % create interpolator
    quadIp = [-(alf - 1) / (alf + 1), 1, (alf - 1) / (alf + 1)];
    Mp = ceil(addPointsAt);
    Mq = floor(N-addPointsAt);

    createTromboneMatrices;
    
    if ~onlyAnalysis
        %% Calculate interpolated pressures
        upMp1 = up(end) * quadIp(3) + wp(1) * quadIp(2) + wp(2) * quadIp(1);
        wpm1 = up(end-1) * quadIp(1) + up(end) * quadIp(2) + wp(1) * quadIp(3);

        %% Calculate v^{n+1/2}
        uvNext(1:end-1) = uv(1:end-1) - lambda / (rho * c) * (up(2:end) - up(1:end-1));
        uvNext(end) = uv(end) - lambda / (rho * c) * (upMp1 - up(end));

        wvNext(2:end) = wv(2:end) - lambda / (rho * c) * (wp(2:end) - wp(1:end-1));
        wvNext(1) = wv(1) - lambda / (rho * c) * (wp(1) - wpm1);

        vNext = v + B * p;

        %% Calculate pressure
        upNext(upRange) = up(upRange) - rho * c * lambda ./ SBar(upRange) .* (SHalf(upRange) .* uvNext(upRange) - SHalf(upRange-1) .* uvNext(upRange-1));
        upNext(1) = up(1) - rho * c * lambda / SBar(1) .* (2 * SHalf(1) * uvNext(1));

        wpNext(wpRange) = wp(wpRange) - rho * c * lambda ./ SBar(wpRange + length(up) - 1) .* (SHalf(wpRange + length(up) - 1) .* wvNext(wpRange+1) - SHalf(wpRange + length(up) - 2) .* wvNext(wpRange));

        if radiation
            wpNext(end) = ((1 - rho * c * lambda * z3) * wp(end) - 2 * rho * c * lambda * (v1 + z4 * p1 - (SHalf(end) .* wvNext(end))/SBar(end))) / (1 + rho * c * lambda * z3);
        end
        v1Next = v1 + k / (2 * Lr) * (wpNext(end) + wp(end));
        p1Next = z1 / 2 * (wpNext(end) + wp(end)) + z2 * p1;

        pNext = radP * p + D * (v + B * p) + radVec;    
    %     pNextTest = [up; wp] + D * vNextTest;    

        %% One step
        pvNextTest = Q * [pv; 1];
        pvNext = pvNextTest(1:end-1);
        if shouldDispCorr && connectedWithP
            displacementCorrection;
        end

        %% Set output from output position
        out(n) = wp(end-1);
    else
        [f, sigma] = analyseQ(full(Q), 1/44100);
        fSave(n, 1:size(f)) = f;
        sigmaSave(n, 1:size(f)) = sigma;

    end
    %% Draw things
    if drawThings && mod (n, drawSpeed) == 0 && n > drawStart
        if onlyAnalysis
            analyseQ(full(Q), 1/44100, false, true);
            drawnow;
        else
            subplot(411)
            hold off;
            plot([uvNext; wvNext])
            hold on;
            plot(vNext);
            plot(pvNext(totSize+1:end));

            subplot(412)
            hold off;
            plot([upNext; wpNext])
            hold on;
            plot(pNext)
            plot(pvNext(1:totSize));

            subplot(413)
            plot(vNext - pvNext(totSize+1:end));

            subplot(414)
            plot(pNext - pvNext(1:totSize));

            if zoomPlot
                for sp = 1:4
                    subplot(4,1,sp)
                    xlim([Mp-10, Mp+10])
                end
            end

            drawnow;
        end
%         pause(0.5)
    end
    
    if ~onlyAnalysis
        %% Update states
        pv = pvNext;

        vPrev = v;
        v = vNext;

        pPrev = p;
        p = pNext;

        uvPrev = uv;
        uv = uvNext;

        wvPrev = wv;
        wv = wvNext;

        upPrev = up;
        up = upNext;

        wpPrev = wp;
        wp = wpNext;
        if radiation
            p1 = p1Next;
            v1 = v1Next;
        end
    end
end
if onlyAnalysis
    plotOneStepTromboneAnalysis;
else
    plot(out)
end