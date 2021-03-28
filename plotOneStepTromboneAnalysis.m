%to be used with modalAnalysis.m
fSave(fSave==0) = nan;
if (Nend-Nstart) < 0
    fSaveRange = 2:size(fSave, 1);
    loopStartRange1 = 2:ceil(numLoops);
else
    fSaveRange = 1:size(fSave,1);
    loopStartRange1 = 1:ceil(numLoops);

end
hold on;
h = plot(real(fSave(fSaveRange, :)));

colours = [];
for colLoop = 1:length(h)
    if mod(colLoop,2) == 0
        colours = [colours; 0,0,1];
    else
        colours = [colours; 1,0,0];
    end
end
%     figure
set(h, {'color', 'Linewidth'}, [num2cell(colours, 2), num2cell(2 * ones(length(h), 1))])
% title ("Modal Analysis $N = " + loopNStart + " \rightarrow" + loopNend + "$", 'interpreter', 'latex');
xlabelsave = num2cell(Nstart:sign(Nend-Nstart):Nend);
set(gca, 'Linewidth', 2, 'Fontsize', 16, 'XTick', loopStartRange1, 'xticklabel', xlabelsave, 'TickLabelInterpreter', 'latex')
xlabel("$N$", 'interpreter', 'latex')
ylabel("Frequency (Hz)", 'interpreter', 'latex')
ylim([0, fs / 2])
