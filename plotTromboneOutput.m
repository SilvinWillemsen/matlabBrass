clc;
clear all;
close all
loadFiles = true;
if loadFiles
%     clear all;
    mode = "Release";
    loadTromboneFiles;
end
lengthSound = length(M);

for n = 1:5:lengthSound
%     subplot(3,1,1)
    hold off;
    plot(1:M(n)+1, pState(n, 1:M(n)+1));
    hold on;    
    plot(M(n)+1:(M(n)+Mw(n) + 1), pState(n, (M(n)+2):(M(n)+Mw(n)+2)));
    pState(n, end)
%     subplot(3,1,2)
% %     plot(energy(1:n, 1:end-1));
% %     legend(["kin", "pot", "rad", "radDamp", "lip", "lipCol", "lipPow", "lipDamp"])
%     subplot(3,1,3)
% % % hold off;
%     plot((sum(energy(1:n, 1:end-1), 2) - energy(n, end))/ 2^(floor(log2(energy(n, end)))));
% % %     hold on
%     plot(scaledTotEnergy(1:n))
% ylim([-3e-16, 3e-16])
    drawnow;
end