clc;
% clear all;
close all
loadFiles = true;
if loadFiles
    clear all;
    onlyLoadFiles = true; % so no plotting when true
    onlyLoadOutput = false;
    drawStart = 1;
    drawSpeed = 100;
    mode = "Debug";
    loadTromboneFiles;
end
lengthSound  = length(output);

% plot(output)
% pause(1)
if onlyLoadFiles
    return;
end
for n = drawStart:drawSpeed:lengthSound
%     subplot(3,1,1)
    if ~onlyLoadOutput
        subplot(2,1,1)
        hold off;
        plot(1:M(n)+1, pState(n, 1:M(n)+1), 'Marker', '.', 'MarkerSize', 10);
        hold on;    
        plot((M(n)+1:(M(n)+Mw(n) + 1)) + alfSave(n), pState(n, (maxM+2):(maxM+Mw(n)+2)), 'Marker', 'o', 'MarkerSize', 2);
%         xlim([M(n) - 5, M(n) + 5])

        subplot(2,1,2)
        hold off;
        plot(1:M(n), vState(n, 1:M(n)), 'Marker', '.', 'MarkerSize', 10);
        hold on;
        plot((M(n)+1:(M(n)+Mw(n))) + alfSave(n), vState(n, (maxM+1):(maxM+Mw(n))), 'Marker', 'o', 'MarkerSize', 2);
        xlim([M(n) - 5, M(n) + 5])
%         subplot(3,1,3)
%         plot(scaledTotEnergy(2:n))
%         plot((sum(energy(2:n, 1:end-1), 2) - energy(2, end))/ 2^(floor(log2(energy(2, end)))));


    else
    %     subplot(2,1,2)
        plot(output(1:n))
    end
%     subplot(3,1,2)
% %     plot(energy(1:n, 1:end-1));
% %     legend(["kin", "pot", "rad", "radDamp", "lip", "lipCol", "lipPow", "lipDamp"])
%     subplot(3,1,3)
% % % hold off;
%     plot((sum(energy(1:n, 1:end-1), 2) - energy(n, end))/ 2^(floor(log2(energy(n, end)))));
% % %     hold on
%     plot(scaledTotEnergy(1:n))
% ylim([-3e-16, 3e-16])
    pause(0.5)
    drawnow;
end
