%% Retrieve new parameters
if changeL
    Ndiff = 1/50;
    Linc = Ndiff * h;
    if (L < Lend)
        L =  L + Linc;
    elseif (L > Lend)
        L = L - Linc;
    end
    
    if (abs(L - Lend) < Linc)
        L = Lend;
    end
    
else
    L = L;
end

LSave(n) = L;

%% Recalculate gridspacing, points, lambda and alpha from new wave speed
% save previous state for comparison later
NPrev = N;

h = c * k; % always stays the same no?
Ninit = L/h;
N = floor(L/h);
Nsave(n) = N;
hSave(n) = h;

lambda = lambdaFact * c * k / h;

alf = (Ninit - N);
alfSave(n) = alf;
% connected with pressures or velocities?
if alternatePV
    connectedWithP = (mod(N,2) == 0);
else
    connectedWithP = true;
end
%     connectedWithP = false;
