clear all;
close all;
drawSpeed = 10;

fs = 44100;
k = 1 / fs;
lengthSound = fs * 2;

f0 = 220;                   % fundamental freq lips
mu = 0.001;                  % mass lips
omega0 = 2 * pi * f0;  % angular freq

%% viscothermal effects
T = 26.85;
[c, rho, eta, nu, gamma] = calcThermoDynConstants(T);
  
sig = 1;
H0 = 0.000;

y = 1;
yPrev = y;

w = 1e-3;
Sr = 1e-3 * w;

S = 0.005;
vNext = 0;
SBar = 0.005;
h = 343 * k;

Pm = 0;
p = 0;
for n = 1:lengthSound
    if y > 0
        theta = 0;
    else
        theta = 1;
    end

    a1 = 2 / k + sig + k * omega0^2;
    a2 = Sr / mu;
    a3 = 2/k * 1/k * (y - yPrev) - omega0^2 * yPrev;
    b1 = S(1) * vNext(1) + h * SBar(1) / (rho * c^2 * k) * (Pm  - p);
    b2 = h * SBar(1) / (rho * c^2 * k);
    c1 = w * subplus(y + H0) * sqrt(2 / rho);
    c2 = b2 + Sr * a2 / a1;
    c3 = b1 - Sr * a3 / a1;
    
    deltaP = sign(c3) * ((-c1 + sqrt(c1^2 + 4 * c2 * abs(c3)))/ (2 * c2))^2;
    
    alpha = 4 / (2 + k * sig + k^2 * omega0^2);
    beta = (k * sig - 2 - k^2 * omega0^2) / (2 + k * sig + k^2 * omega0^2);
    epsilon = 2 * k^2 * Sr / (mu * (2 + k*sig + k^2 * omega0^2));
    
    yNext(n) = alpha * y + beta * yPrev;% + epsilon * deltaP;
%     yNext(n) = ((2 * mu / k^2 - (1 + 3 * theta) * K) * y + ((1 + 4 * theta) * R / (2 * k) - mu / k^2) * yPrev) / (mu / k^2 + (1 + 4 * theta) * R / (2 * k));
    Ub = w * subplus(y + H0) * sign(deltaP) * sqrt(2 * abs(deltaP)/rho);
    Ur = Sr * 1/(2*k) * (yNext(n) - yPrev);
%     zNext(n) = (2 * z - zPrev - (1 + 3 * theta) * omegaSq * k^2 * z + (1 + 4 * theta) * R * k^2 / 2 * zPrev) / (1 + (1 + 4 * theta) * R * k^2 / 2);
    hReed(n) = mu / 2 * ((1/k * (y - yPrev))^2 + omega0^2 * (y^2 + yPrev^2) / 2);
    rOChReed(n) = mu * (1/(2*k) * (yNext(n) - yPrev) * 1/k^2 * (yNext(n) - 2 * y + yPrev) ...
         + omega0^2 * 1/(2*k) * (yNext(n) - yPrev) * 1/2 * (yNext(n) + yPrev));
    qReed(n) = mu * sig * (1/(2*k) * (yNext(n) - yPrev))^2;% + w * subplus(y + H0) * (sqrt(2 / rho) * abs(deltaP)^(3/2));
    qHReed(n) = k * sum (qReed);
    pReed(n) = -(Ub + Ur) * Pm;
    
    totEnergy(n) = hReed(n) + k * sum(qReed);
    totROCEnergy(n) = rOChReed(n) + qReed(n);
    if mod(n, drawSpeed) == 0
        subplot(311)
        plot(yNext(1:n))
        subplot(312)
        hold off;
        plot(hReed(1:n) - hReed(1));
        hold on;
        plot(qHReed(1:n))
%         plot(qReed(1:n))
        subplot(313)
        plot(totEnergy(1:n))
        drawnow;
    end
    
    yPrev = y;
    y = yNext(n);
end