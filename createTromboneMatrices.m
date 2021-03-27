totSize = Mp+1+Mq+1;

Bu = sparse(1:Mp+1, 1:Mp+1, ones(1, Mp+1) * lambda / (rho * c), Mp+1, Mp+1) + ...
     sparse(1:Mp, 2:Mp+1, ones(1, Mp) * -lambda / (rho * c), Mp+1, Mp+1);
 
Bw = sparse(1:Mq+1, 1:Mq+1, ones(1, Mq+1) * -lambda / (rho * c), Mq+1, Mq+1) + ...
     sparse(2:Mq+1, 1:Mq, ones(1, Mq) * lambda / (rho * c), Mq+1, Mq+1);

B = zeros(totSize);
B(1:Mp+1, 1:Mp+1) = Bu;
B(Mp+2:end, Mp+2:end) = Bw;

%% add interpolation
B(Mp+1, Mp+1:Mp+3) = B(Mp+1, Mp+1:Mp+3) - lambda / (rho * c) * fliplr(quadIp);
B(Mp+2, Mp:Mp+2) = B(Mp+2, Mp:Mp+2) + lambda / (rho * c) * quadIp;
% B(M+1, M+3) = -lambda / (rho * c);
% B(M+2, M) = lambda / (rho * c);
%% matices to calculate pNext
Du = sparse(1:Mp+1, 1:Mp+1, ones(1, Mp+1) .* -rho .* c .* lambda .* SHalf(1:Mp+1)'./SBar(1:Mp+1)', Mp+1, Mp+1) + ...
     sparse(2:Mp+1, 1:Mp, ones(1, Mp) .* rho .* c .* lambda .* SHalf(1:Mp)'./SBar(2:Mp+1)', Mp+1, Mp+1);
 
Du(1, 1) = 2 * Du(1, 1);

Dw = sparse(1:Mq+1, 1:Mq+1, ones(1, Mq+1) .* rho .* c .* lambda .* SHalf(Mp:end)'./SBar(Mp+1:end)', Mq+1, Mq+1) + ...
     sparse(1:Mq, 2:Mq+1, ones(1, Mq) .* -rho .* c .* lambda .* SHalf(Mp+1:end)'./SBar(Mp+1:end-1)', Mq+1, Mq+1);

radMat = zeros(totSize, 1);
radP = sparse(1:totSize, 1:totSize, ones(1, totSize), totSize, totSize);
if radiation
    radP(end) = (1 - rho * c * lambda * z3) / (1 + rho * c * lambda * z3);
    radMat(end) = - 2 * rho * c * lambda * (v1 + z4 * p1) / (1 + rho * c * lambda * z3);
    Dw(end, end) = (2 * rho * c * lambda * SHalf(end)/SBar(end)) / (1 + rho * c * lambda * z3);
    
%     wpNext(end) = ((1 - rho * c * lambda * z3) * wp(end) - 2 * rho * c * lambda * (v1 + z4 * p1 - SHalf(end)/SBar(end) .* wvNext(end))) / (1 + rho * c * lambda * z3);
else
    Dw(end, end) = 0;
end

D = zeros(totSize);
D(1:Mp+1, 1:Mp+1) = Du;
D(Mp+2:end, Mp+2:end) = Dw;

%% onestep form

% A1 = [eye(totSize), D; zeros(totSize), eye(totSize)];
% B1 = [eye(totSize), zeros(totSize); B, eye(totSize)];
Q = [radP + D * B,             D,                sparse(1:totSize, 1:totSize, radMat', totSize, totSize); 
     B,                        eye(totSize),     zeros(totSize);
     eye(totSize),             zeros(totSize),   zeros(totSize)];

