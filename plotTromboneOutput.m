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

for n = 1:2:lengthSound
    subplot(2,1,1)
    hold off;
    plot(1:M(n), pState(n, 1:M(n)));
    hold on;
    plot(M(n):(M(n)+Mw(n) + 1), pState(n, (M(n)+1):end));
    pState(n, end)
    subplot(2,1,2)
    plot(energy(1:n));
ylim([-3e-16, 3e-16])
    drawnow;
end