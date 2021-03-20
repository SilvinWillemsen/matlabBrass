%{
    Full trombone: lip, dynamic tube with proper geometry, radiation 
%}

% clear all;
close all;

fs = 44100;             % Sample rate (Hz)
k = 1/fs;               % Time step (s)
lengthSound = fs;       % Duration (s)

% drawing variables
drawThings = false;
drawSetting = 0; % change to 0 to go back to previous drawing
zoomPlot = false;

drawSpeed = 1000;
drawStart = 0;
drawSpeedInit = drawSpeed;

fixedNonInterpolatedL = false;

changeL = ~fixedNonInterpolatedL;
radiation = true;

connectedToLip = true;

LnonExtended = 2.658;
Lextended = 3.718;

%% viscothermal effects
T = 26.85;
[c, rho, eta, nu, gamma] = calcThermoDynConstants (T);

%% Tube variables
h = c * k;              % Grid spacing (m)
NnonExtended = floor(LnonExtended / h);
Nextended = floor(Lextended / h);

Nstart = NnonExtended;
Nend = Nextended;

if fixedNonInterpolatedL
    L = floor(Nstart) * h;
    Ninit = L / h;
else
    Ninit = Nstart;
end
lambdaFact = 0.99999;
lambda = lambdaFact * c * k / h      % courant number

LInit = Ninit*h;
Lend = Nend * h;

L = LInit;
Ninit = L/h;
N = floor(Ninit);         % Number of points (-)
alf = Ninit - N;

%% Lip Collision
Kcol = 10000;
alfCol = 3; 

%% Set cross-sectional geometry
setToOnes = false;
[S, SHalf, SBar, addPointsAt] = setTube (N+1, NnonExtended, 0, setToOnes);

% Quick note: N is the number of spaces between the points so the number of points is N+1

%% Lip variables
f0 = 300;

Mlip = 5.37e-5;                % mass lips
omega0 = 2 * pi * f0;   % angular freq

sig = 5;                % damping
H0 = 2.9e-4;                % equilibrium
y = 0;                      % initial lip state
yNext = zeros(lengthSound, 1);                  
deltaP = 0;
psi = 0;
Ub = 0;
Ur = 0;
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
upPrev = zeros(ceil(addPointsAt) + 1, 1);
wpNext = zeros(floor(N-addPointsAt) + 1, 1);
wp = zeros(floor(N-addPointsAt) + 1, 1);
wpPrev = zeros(floor(N-addPointsAt) + 1, 1);

connectedWithP = (mod(N,2) == 0);
% connectedWithP = false;

if connectedWithP
    uvNext = zeros(ceil(addPointsAt), 1); 
    uv = zeros(ceil(addPointsAt), 1);
    wvNext = zeros(floor(N-addPointsAt), 1); % the total v vector is one smaller than the total p vector
    wv = zeros(floor(N-addPointsAt), 1);

else
    uvNext = zeros(ceil(addPointsAt) + 1, 1); 
    uv = zeros(ceil(addPointsAt) + 1, 1);
    wvNext = zeros(floor(N-addPointsAt) + 1, 1); % the total v vector is one smaller than the total p vector
    wv = zeros(floor(N-addPointsAt) + 1, 1);

end
         
uvMph = 0;
uvNextMph = 0;
wvmh = 0;
wvNextmh = 0;

upMp1Prev = 0;
wpm1Prev = 0;

if ~connectedToLip
    exciteU = false;
    if exciteU
        width = floor(0.2 * length(up));
        loc = floor(0.5 * length(up));
        inputRange = (loc-width) : (loc+width);
        up(inputRange) = up(inputRange) + hann(width*2+1);
        upPrev = up;
    else
        width = floor(0.2 * length(wp));
        loc = floor(0.5 * length(wp));
        inputRange = (loc-width) : (loc+width);
        wp(inputRange) = wp(inputRange) + hann(width*2+1);
        wpPrev = wp;
    end
        
end

% Initialise output
out = zeros (lengthSound, 1);

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
p1Prev = 0;
v1Prev = 0;

Pmprev = 0;
flag = false;
statesSave = [];

for n = 1:lengthSound
    
    retrieveAndRecalculateParams;
    addAndRemovePoints;
    
    % create interpolator
    quadIp = [-(alf - 1) / (alf + 1), 1, (alf - 1) / (alf + 1)];
    
    if connectedWithP
        %% Calculate interpolated pressures
        upMp1 = up(end) * quadIp(3) + wp(1) * quadIp(2) + wp(2) * quadIp(1);
        wpm1 = up(end-1) * quadIp(1) + up(end) * quadIp(2) + wp(1) * quadIp(3);

        %% Calculate v^{n+1/2}
        uvNext = uv - lambda / (rho * c) * (up(2:end) - up(1:end-1));
        uvNextMph = uvMph - lambda / (rho * c) * (upMp1 - up(end));

        wvNext = wv - lambda / (rho * c) * (wp(2:end) - wp(1:end-1));
        wvNextmh = wvmh - lambda / (rho * c) * (wp(1) - wpm1);
        
%         wvNextmh - uvNext(end)
%         uvNextMph - wvNext(1)
%         up(end) - wp(1)

    else
        %% Calculate interpolated velocities
        uvMph = uv(end) * quadIp(3) + wv(1) * quadIp(2) + wv(2) * quadIp(1);
        wvmh = uv(end-1) * quadIp(1) + uv(end) * quadIp(2) + wv(1) * quadIp(3);

        %% Calculate p^n
        up(upRange) = upPrev(upRange) - rho * c * lambda ./ SBar(upRange) .* (SHalf(upRange) .* uv(upRange) - SHalf(upRange-1) .* uv(upRange-1));
        upMp1 = upMp1Prev - rho * c * lambda / SBar(length(up)) * (SHalf(length(up)) .* uvMph - SHalf(length(up) - 1) .* uv(end));
%         if connectedToLip
            up(1) = upPrev(1) - rho * c * lambda / SBar(1) .* (-2 * (Ub + Ur) + 2 * SHalf(1) * uv(1));
%         end

        wp(wpRange-1) = wpPrev(wpRange-1) - rho * c * lambda ./ SBar(wpRange + length(up) - 2) .* (SHalf(wpRange + length(up) - 2) .* wv(wpRange) - SHalf(wpRange + length(up) - 3) .* wv(wpRange-1));
        wpm1 = wpm1Prev - rho * c * lambda ./ SBar(length(up)) .* (SHalf(length(up)) .* wv(1) - SHalf(length(up) - 1) .* wvmh);
        
        if radiation
            wp(end) = ((1 - rho * c * lambda * z3) * wpPrev(end) - 2 * rho * c * lambda * (v1 + z4 * p1 - (SHalf(end) .* wv(end))/SBar(end))) / (1 + rho * c * lambda * z3);
        end
        
        v1Next = v1 + k / (2 * Lr) * (wp(end) + wpPrev(end));
        p1Next = z1 / 2 * (wp(end) + wpPrev(end)) + z2 * p1;

        %% Calculate v^{n+1/2}
        uvNext(1:end-1) = uv(1:end-1) - lambda / (rho * c) * (up(2:end) - up(1:end-1));
        uvNext(end) = uv(end) - lambda / (rho * c) * (upMp1 - up(end));
        
        wvNext(2:end) = wv(2:end) - lambda / (rho * c) * (wp(2:end) - wp(1:end-1));
        wvNext(1) = wv(1) - lambda / (rho * c) * (wp(1) - wpm1);
        
%         wpm1 - up(end)
%         upMp1 - wp(1)


    end
    
    Pm = amp;
    
    if connectedToLip
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
    else
        Ub = 0;
        Ur = 0;
    end
    
%     retrieveAndRecalculateParams;
%     
%     addAndRemovePoints;    
    
    if connectedWithP
        %% Calculate pressure
        upNext(upRange) = up(upRange) - rho * c * lambda ./ SBar(upRange) .* (SHalf(upRange) .* uvNext(upRange) - SHalf(upRange-1) .* uvNext(upRange-1));
        upNext(1) = up(1) - rho * c * lambda / SBar(1) .* (-2 * (Ub + Ur) + 2 * SHalf(1) * uvNext(1));
        upNext(end) = up(end) - rho * c * lambda ./ SBar(length(up)) .* (SHalf(length(up)) .* uvNextMph - SHalf(length(up) - 1) .* uvNext(end));

        wpNext(wpRange) = wp(wpRange) - rho * c * lambda ./ SBar(wpRange + length(up) - 2) .* (SHalf(wpRange + length(up) - 2) .* wvNext(wpRange) - SHalf(wpRange + length(up) - 3) .* wvNext(wpRange-1));
        wpNext(1) = wp(1) - rho * c * lambda ./ SBar(length(up)) .* (SHalf(length(up)) .* wvNext(1) - SHalf(length(up) - 1) .* wvNextmh);
        if radiation
            wpNext(end) = ((1 - rho * c * lambda * z3) * wp(end) - 2 * rho * c * lambda * (v1 + z4 * p1 - (SHalf(end) .* wvNext(end))/SBar(end))) / (1 + rho * c * lambda * z3);
        end
        v1Next = v1 + k / (2 * Lr) * (wpNext(end) + wp(end));
        p1Next = z1 / 2 * (wpNext(end) + wp(end)) + z2 * p1;
    end
    %% Set output from output position
    out(n) = wp(end-1);

    %% Draw things
    if drawThings && mod (n, drawSpeed) == 0 && n > drawStart

        locsLeft = 0:length(up)-1;
        locsRight = (0:length(wp)-1)+length(up) + alf; 
        if connectedWithP
            locsRight = locsRight - 1;
        end
        plotPrev = false;
        if plotPrev
              %% Plot pressures
            subplot(2,1,1)
            if connectedWithP
                hold off;
                plot(locsLeft, upPrev, '-o');
                hold on;
                plot(locsRight, wpPrev, '-o');
            else
                hold off;
                plot([locsLeft, locsLeft(end) + 1], [upPrev; upMp1Prev], '-o');
                hold on;
                plot([locsRight(1) - 1, locsRight], [wpm1Prev; wpPrev], '-o');
                plot((locsLeft(end) + 1), upMp1Prev, 'Marker', 'o', 'MarkerSize', 10, 'Color', 'r');
                plot((locsRight(1) - 1), wpm1Prev, 'Marker', 'o', 'MarkerSize', 10,  'Color', 'b');

            end
    %         xlim([locsLeft(end-10), locsRight(10)])
%             ylim([-2, 2])
            %% Plot velocities
            subplot(2,1,2)
            if connectedWithP
                hold off;
                plot(locsLeft + 0.5, [uv; uvMph], 'Marker', '.', 'MarkerSize', 10, 'Color', 'r');
                hold on;
                plot(locsRight - 0.5, [wvmh; wv], 'Marker', '.', 'MarkerSize', 10,  'Color', 'b');
                plot(locsLeft(end) + 0.5, uvMph, 'Marker', 'o', 'MarkerSize', 10, 'Color', 'r');
                plot(locsRight(1) - 0.5, wvmh, 'Marker', 'o', 'MarkerSize', 10,  'Color', 'b');
            else
                hold off;
                plot(locsLeft + 0.5, uv, 'Marker', '.', 'MarkerSize', 10, 'Color', 'r');
                hold on;
                plot(locsRight - 0.5, wv, 'Marker', '.', 'MarkerSize', 10,  'Color', 'b');

            end
%             ylim([-1e-2, 1e-2])

    %         xlim([locsLeft(end-10), locsRight(10)])
            pause(0.1)
            drawnow;
        else

            %% Plot pressures
            subplot(2,1,1)
            if connectedWithP
                hold off;
                plot(locsLeft, up, '-o');
                hold on;
                plot(locsRight, wp, '-o');
            else
                hold off;
                plot([locsLeft, locsLeft(end) + 1], [up; upMp1], '-o');
                hold on;
                plot([locsRight(1) - 1, locsRight], [wpm1; wp], '-o');
                plot((locsLeft(end) + 1), upMp1, 'Marker', 'o', 'MarkerSize', 10, 'Color', 'r');
                plot((locsRight(1) - 1), wpm1, 'Marker', 'o', 'MarkerSize', 10,  'Color', 'b');

            end
            if zoomPlot
                xlim([locsLeft(end-10), locsRight(10)])
%                     xlim([locsRight(end-10), locsRight(end)])

            end
%             ylim([-2, 2])
            %% Plot velocities
            subplot(2,1,2)
            if connectedWithP
                hold off;
                plot(locsLeft + 0.5, [uvNext; uvNextMph], 'Marker', '.', 'MarkerSize', 10, 'Color', 'r');
                hold on;
                plot(locsRight - 0.5, [wvNextmh; wvNext], 'Marker', '.', 'MarkerSize', 10,  'Color', 'b');
                plot(locsLeft(end) + 0.5, uvNextMph, 'Marker', 'o', 'MarkerSize', 10, 'Color', 'r');
                plot(locsRight(1) - 0.5, wvNextmh, 'Marker', 'o', 'MarkerSize', 10,  'Color', 'b');
            else
                hold off;
                plot(locsLeft + 0.5, uvNext, 'Marker', '.', 'MarkerSize', 10, 'Color', 'r');
                hold on;
                plot(locsRight - 0.5, wvNext, 'Marker', '.', 'MarkerSize', 10,  'Color', 'b');

            end
%             ylim([-1e-2, 1e-2])
            if zoomPlot
                xlim([locsLeft(end-10), locsRight(10)])
%                     xlim([locsRight(end-10), locsRight(end)])
            end
            pause(0.5)
            drawnow;
        end
    end

    %% Update states
    uv = uvNext;
    wv = wvNext;
    
    upPrev = up;
    if connectedWithP
        up = upNext;
    end
    
    wpPrev = wp;     
    if connectedWithP
        wp = wpNext;
    end
    
    upMp1Prev = upMp1;
    wpm1Prev = wpm1;
    
    uvMph = uvNextMph;
    
    wvmh = wvNextmh;
    
    if radiation
        p1 = p1Next;
        v1 = v1Next;
    end
    
    yPrev = y;
    y = yNext(n);
    psiPrev = psi; 

end
plot(out)
