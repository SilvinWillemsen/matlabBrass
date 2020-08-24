clear all;
close all;
drawSpeed = 10000;

fs = 44100;
k = 1 / fs;
lengthSound = fs * 10;

f0 = 220;                   % fundamental freq lips
m = 0.001;                  % mass lips
omegaSq = (2 * pi * f0)^2;  % angular freq

%% viscothermal effects
tempAbove2685 = 0;
c = 347.23 * (1 + 0.00166 * tempAbove2685);         % speed of sound in air (m/s)
rho = 1.1769 * (1 - 0.00335 * tempAbove2685);       % air density (kg / m^3)
eta = 1.846 * 1e-5 * (1 + 0.0025 * tempAbove2685);  % shear viscosity
nu = 0.8410 * (1 - 0.0002 * tempAbove2685);         % root of prandtl number
gamma = 1.4017 * (1 - 0.00002 * tempAbove2685);     % ratio of specific heats
    
K = omegaSq * m;            % spring constant lip
Pm = 20;                    % Blowing pressure (kPa)
Rspring = 0.23;             % initial damping
x0 = (0.017 * Pm / 100 + 0.23) * 1e-3;       % mean displacement
q = 2.0;                    % exponent value
Ab = 1e-6;                  % area bottom (m^2)

wR = 0.33;
z = 0;
zPrev = z;
prevP = 1;
for n = 1:lengthSound
    if z > 0
        theta = 0;
    else
        theta = 1;
    end
%     if n > lengthSound / 2
%         Pm = 0;
%     end

    R = max (Rspring - wR * abs (x0 / z)^q, 0);
    Rsave(n) = R;
    
    A = c * (1 - theta) * sqrt(2 * rho);
    p1 = Pm + ((-A + sqrt(A^2 - 4 * Pm)) / 2)^2;
    p2 = Pm - ((-A + sqrt(A^2 + 4 * Pm)) / 2)^2;
    
    if sign(p1) == 1 && sign(p2) == 0
        p = p1;
        prevP = 1;
    elseif sign(p2) == 1 && sign(p1) == 0
        p = p2;
        prevP = 2;
    else % if both are positive
        if prevP == 1
            p = p1;
        else
            p = p2;
        end
    end
    pSave(n) = p;
    Fsave(n) = gamma * (Pm - p) + Ab * (1 - theta) * p;
    zNext(n) = (Fsave(n) + (2 * m / k^2 - (1 + 3 * theta) * K) * z + ((1 + 4 * theta) * R / (2 * k) - m / k^2) * zPrev) / (m / k^2 + (1 + 4 * theta) * R / (2 * k));

%     zNext(n) = (2 * z - zPrev - (1 + 3 * theta) * omegaSq * k^2 * z + (1 + 4 * theta) * R * k^2 / 2 * zPrev) / (1 + (1 + 4 * theta) * R * k^2 / 2);
    
    if mod(n, drawSpeed) == 0
        plot(zNext(1:n))
        drawnow;
    end
    
    zPrev = z;
    z = zNext(n);
end