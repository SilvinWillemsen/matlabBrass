clear all;
close all;

fs = 44100;
lengthSound = fs;

k = 1/fs;
c = 300;

h = c*k;
N = floor(1/h);
h = 1/N;

lambdaSq = (c * k / h)^2;

sig = 100;

uNext = zeros(N, 1);
u = zeros(N, 1);
u(floor(N/2-5):floor(N/2+5)) = hann(11);
uPrev = u;

range = 3:N-2;
potEnergyRange = 2:N-1;
dampEnergy = zeros(lengthSound, 1);
for n = 1:lengthSound
    uNext(range) = (2 * u(range) - uPrev(range) + lambdaSq * (u(range+1) - 2 * u(range) + u(range-1)) + sig * k * uPrev(range)) / (1 + sig*k);
    
    
    kinEnergy(n) = 1/2 * h * sum((1/k * (u - uPrev)).^2);
    potEnergy(n) = c^2 / 2 * h * 1/h^2 * sum((u(potEnergyRange+1) - u(potEnergyRange)) .* (uPrev(potEnergyRange+1) - uPrev(potEnergyRange)));
    rOCdampEnergy(n) = 2 * sig * h * sum ((1/(2*k) * (uNext - uPrev)).^2);
    
    idx = n - (1 * (n~=1));
    dampEnergy(n) = dampEnergy(idx) + k * rOCdampEnergy(n);
    totEnergy(n) = kinEnergy(n) + potEnergy(n) + dampEnergy(idx);
    
    %drawthings
    subplot(211)
    plot(u)
    
    subplot(212)
    if n>10
        plot(totEnergy(10:n) / totEnergy(10) - 1)
    end
    drawnow;

    
    uPrev = u;
    u = uNext;
end