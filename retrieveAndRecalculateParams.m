%% Retrieve new parameters
if changeL
    Linc = 0.00002;
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

% connected with pressures or velocities?
connectedWithP = (mod(N,2) == 0);
%     connectedWithP = false;
