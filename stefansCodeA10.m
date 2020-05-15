close all;
clear all;

% matlab script vocalfd.m
% finite difference vocal tract simulation % radiation loss included
% simple glottal source waveform
% static pitch

%%%%% real-time plotting variables

drawThings = true; % disable if you don't want a real-time plot
drawSpeed = 3; 

%%%%%% begin global parameters

SR = 44100;         % sample rate (Hz)
L = 1;              % tract length (m)
S0 = 0.00025;       % vocal tract surface area, left end (m^2)
c = 340;            % wave speed (m/s)
f0 = 120;           % fundamental frequency (Hz)
TF = 1;             % simulation duration (s)
% vocal tract profile, non-dimensional [pos S] pairs
    % /E/
        S = [0 1;0.09 0.4;0.11 2.4;0.24 2.4;0.26 3.2;0.29 3.2;0.32 4.2; ... 
            0.41 4.2;0.47 3.2;0.59 1.8;0.65 1.6;0.71 1.6;0.74 1;0.76 0.8;... 
            0.82 0.8;0.88 2;0.91 2;0.94 3.2;1 3.2];
    % /A/
    % S = [0 1;0.03 0.60;0.09 0.4;0.12 1.6;0.18 0.6;0.29 0.2;0.35 0.4;...
    % 0.41 0.8;0.47 1;0.50 0.6;0.59 2;0.65 3.2;0.85 3.2;0.94 2;1 2];
    
%%%%% end global parameters
    
% begin derived parameters

k = 1/SR;           % time step
NF = floor(TF*SR);  % sample duration
gamma = c/L;        % stability condition/scheme parameters

h = gamma*k; N = floor(1/h); h = 1/N; lambda = gamma*k/h;
S = interp1(S(:,1),S(:,2),[0:h:1])';    % interpolate vocal tract profil
alf = 2.0881*L*sqrt(1/(S0*S(N+1)));     % radiation parameter
% bet = 0.7407/gamma;                   % radiation parameter (set to zero to not have any losses in the system)
bet = 0;
Sav = [S(1); 0.25*(S(3:N+1)+2*S(2:N)+S(1:N-1)); S(N+1)]; 
Sr = 1.5*S(N+1)-0.5*S(N);
sr = 0.5*lambda^2*((S(2:N)+S(3:N+1))./Sav(2:N)); 
sl = 0.5*lambda^2*((S(2:N)+S(1:N-1))./Sav(2:N)); 
s0 = 2*(1-lambda^2);
q1 = alf*gamma^2*k^2*Sr/(Sav(N+1)*h); q2 = bet*gamma^2*k*Sr/(Sav(N+1)*h); 
r1 = 2*lambda^2/(1+q1+q2); r2 = -(1+q1-q2)/(1+q1+q2);
g1 = -(k^2*gamma^2/h/S(1))*(3*S(1)-S(2));

% initialize grid functions and output

% I removed input and excited with an impulse instead to be better able to look at the energy
Psi = zeros(N+1,1); 
Psi1 = zeros(N+1,1);
Psi1(floor(N/4) : floor(N/4) + 9) = hann(10);
Psi2 = Psi1; 
% uin = sin(2*pi*[0:NF-1]*k*f0);
% uin = 0.5*(uin+abs(uin));
out = zeros(NF,1);

%%%%% set energy scaling
scaling = ones(N+1,1);
epsilonL = 1/2 * (S(1) + S(2)) / Sav(1);
epsilonR = 1/2 * (S(N+1) + S(N)) / Sav(N+1);
scaling(N+1) = epsilonR / 2;
scaling(1) = epsilonL / 2;

%%%%% initialise energies
kinEnergy = zeros(NF, 1);
potEnergy = zeros(NF, 1);
boundaryEnergy = zeros(NF, 1);

%%%%%% begin main loop

for n=1:NF
    % update scheme without boundaries
    Psi(2:N) = s0*Psi1(2:N)+sl.*Psi1(1:N-1)+sr.*Psi1(3:N+1)-Psi2(2:N);
    
    % Left boundary (added a missing s0*Psi1(N+1)/(1+q1+q2) term below) 
    Psi(N+1) = s0*Psi1(N+1)/(1+q1+q2) + r1*Psi1(N) + r2*Psi2 (N+1);
    
    % Right boundary (removed input to better look at the energy plot) 
    Psi(1) = s0*Psi1(1)+2*lambda^2*Psi1(2)-Psi2(1);%+g1*uin(n);
    out(n) = SR*(Psi(N+1)-Psi1(N+1));
    
    %% energies
    kinEnergy(n) = 1/2 * sum(h * Sav .* scaling .* (1/k * (Psi1 - Psi2)).^2);
    potEnergy(n) = -h * gamma^2 / 2 * 1/h^2 * sum((S(1:N) + S(2:N+1)) * 0.5...
        .* (Psi1(2:N+1) - Psi1(1:N)) .* (Psi2(2:N+1) - Psi2(1:N)));
    boundaryEnergy(n) = 2 * (1 - epsilonR / 2) * 0.5 * (S(N+1) + S(N)) * gamma^2 * alf / 4 * (Psi1(N+1)^2 + Psi2(N+1)^2);
    totEnergy(n) = kinEnergy(n) - potEnergy(n) + boundaryEnergy(n);
    
    %%%% draw plots
    if drawThings && mod(n, drawSpeed) == 0
        subplot(2,1,1)
        hold off;
        plot(Psi, 'b', 'Linewidth', 2)
        hold on;
        plot(sqrt(S),'k');
        plot(-sqrt(S),'k');
        title("$\Psi$", 'interpreter', 'latex')
        xlim([1,N+1])
        set(gca, 'Fontsize', 16)

        subplot(2,1,2)
        plot(totEnergy(1:n) / totEnergy(1) - 1)
        title("Normalised total energy")
        set(gca, 'Fontsize', 16)
        drawnow;
    end
    Psi2 = Psi1; Psi1 = Psi;
end
%%%%%% end main loop
% plot vocal tract profile and output spectrum
subplot(2,1,1); 
plot([0:h:1], sqrt(S),'k', [0:h:1], -sqrt(S),'k');
title('Vocal Tract Profile'); 
xlabel('x'); 
ylabel('sqrt(S)'); 
subplot(2,1,2); 
plot([0:NF-1]*SR/NF, 10*log10(abs(fft(out))), 'k'); 
title('Output Spectrum'); xlabel('f');
ylabel('pressure (dB)');
% play sound
soundsc(out, SR);