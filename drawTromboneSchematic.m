%{ 
    Script to draw the Trombone Schematic used for the DAFx21 paper
%}
close all;
clear all;

figure("Position", [0, 1000, 1300, 350])
set(gca, 'Position', [0.00222222222222222 0.01 0.994444444444444 0.986666666666667]);

%% options
dashedLineWidth = 1;
timeIndices = false; 

fontSize = 18;

numPLeft = 4;
numPRight = 3; 

numQRight = 1;
numQLeft = 3; 

alignDashedToOutside = false;
curvatureArcs = 0.8 * pi;
drawVnmh = true;

alf = 0.25;

%% prepare figure locations
vScale = 1; % determines curvature of arrows
xPLocs = 0:(numPLeft - 1);
xPLocs = [xPLocs, ((0.5:0.5:1.5)+xPLocs(end))];
xPLocs = [xPLocs, ((1:numPRight-1)+xPLocs(end))];

xQLocs = 0:(numQLeft - 1);
xQLocs = [xQLocs, ((0.5:0.5:1.5)+xQLocs(end))];
xQLocs = [xQLocs, ((1:numQRight-1)+xQLocs(end))];

xQLocs = xQLocs + xPLocs(end) + alf;

%% draw time indices and the grid lines
timeSteps = ["$n+1$", "$n+1/2$", "$n$", "$n-1/2$"];
j = 1;
if drawVnmh
    endidx = -3;
else
    endidx = -2
end

for nIdx = 0:-1:endidx
    plot([-0.25, xQLocs(end)+0.25], [nIdx, nIdx] * vScale, 'color', [0.5, 0.5, 0.5]);
    text(-0.3, nIdx * vScale, timeSteps(j), 'interpreter', 'latex', 'Fontsize', 16, 'horizontalAlignment', 'right')
    hold on;
    j = j + 1;
end

% set limits
xlim([-0.75, xQLocs(end) + 0.25])
if drawVnmh
    ylim([-3.5, 0.4] * vScale)  
else
    ylim([-2.5, 0.4] * vScale)
end

%% Draw arrows showing dataflow 
arrowLineWidth = 1.5;
arrowHeadWidth = 5;
arrowHeadLength = 10;
includeStraightArrows = false;
alphArr = 0.8;
alphArrInterp = 0.4;
tromboneFigDrawArrows;

for n = 0:1
    %% plot pressures

    % p
    plot([xPLocs(1), xPLocs(numPLeft+1)], [0, 0] - n * 2 * vScale, 'LineWidth', 2, 'Color', 'b');
    plot([xPLocs(numPLeft+1), xPLocs(numPLeft+2)], [0, 0] - n * 2 * vScale, '--', 'LineWidth', dashedLineWidth, 'Color', 'b');
    plot([xPLocs(numPLeft+2), xPLocs(end)], [0, 0] - n * 2 * vScale, 'LineWidth',2, 'Color', 'b');

    pXLocs = [xPLocs(1:numPLeft), xPLocs(end-numPRight+1:end)];
    scatter(pXLocs, zeros(1,length(pXLocs)) - n * 2 * vScale, 400, 'b', '.');

    % q
    plot([xQLocs(1), xQLocs(numQLeft+1)], [0, 0] - n * 2 * vScale, 'LineWidth', 2, 'Color', 'r');
    plot([xQLocs(numQLeft+1), xQLocs(numQLeft+2)], [0, 0] - n * 2 * vScale, '--', 'LineWidth', dashedLineWidth, 'Color', 'r');
    plot([xQLocs(numQLeft+2), xQLocs(end)], [0, 0] - n * 2 * vScale, 'LineWidth', 2, 'Color', 'r');

    qXLocs = [xQLocs(1:numQLeft), xQLocs(end-numQRight+1:end)];

    scatter(qXLocs, zeros(1,length(qXLocs)) - n * 2 * vScale, 400, 'r', '.');

    %% plot velocities
    vOffset = -1;
    if n ~= 1 || drawVnmh

        % v
        if alignDashedToOutside
            xVLocs = xPLocs(2:end) - 0.5;
            plot([xVLocs(1), xVLocs(numPLeft)], [0, 0] + (vOffset - n * 2) * vScale, 'LineWidth', 2, 'Color', 'b');
            plot([xVLocs(numPLeft), xVLocs(numPLeft+1)], [0, 0] + (vOffset - n * 2) * vScale, '--', 'LineWidth', dashedLineWidth, 'Color', 'b');
            plot([xVLocs(numPLeft+1), xVLocs(end)], [0, 0] + (vOffset - n * 2) * vScale, 'LineWidth', 2, 'Color', 'b');

            vXLocs = [xVLocs(1:numPLeft-1), xVLocs(end-numPRight+1:end)];
        else
            xVLocs = xPLocs(1:end-1) + 0.5;
            plot([xVLocs(1), xVLocs(numPLeft+1)], [0, 0] + (vOffset - n * 2) * vScale, 'LineWidth', 2, 'Color', 'b');
            plot([xVLocs(numPLeft+1), xVLocs(numPLeft+2)], [0, 0] + (vOffset - n * 2) * vScale, '--', 'LineWidth', dashedLineWidth, 'Color', 'b');
            plot([xVLocs(numPLeft+2), xVLocs(end)], [0, 0] + (vOffset - n * 2) * vScale, 'LineWidth', 2, 'Color', 'b');

            vXLocs = [xVLocs(1:numPLeft), xVLocs(end-numPRight+2:end)];
        end

        scatter(vXLocs, zeros(1, length(vXLocs)) + (vOffset - n * 2) * vScale,  400, 'b', '.');
        scatter(vXLocs(end) + 1, 0 + (vOffset - n * 2) * vScale, 80, 'b', 's', 'LineWidth', 1.5);


        % w

        if alignDashedToOutside
            xWLocs = xQLocs(1:end-1) + 0.5;
            plot([xWLocs(1), xWLocs(numQLeft+1)], [0, 0] + (vOffset - n * 2) * vScale, 'LineWidth', 2, 'Color', 'r');
            plot([xWLocs(numQLeft+1), xWLocs(numQLeft+2)], [0, 0] + (vOffset - n * 2) * vScale, '--', 'LineWidth', dashedLineWidth, 'Color', 'r');
            plot([xWLocs(numQLeft+2), xWLocs(end)], [0, 0] + (vOffset - n * 2) * vScale, 'LineWidth', 2, 'Color', 'r');

            wXLocs = [xWLocs(1:numQLeft), xWLocs(end-numQRight+2:end)];
        else
            xWLocs = xQLocs(2:end) - 0.5;
            plot([xWLocs(1), xWLocs(numQLeft)], [0, 0] + (vOffset - n * 2) * vScale, 'LineWidth', 2, 'Color', 'r');
            plot([xWLocs(numQLeft), xWLocs(numQLeft+1)], [0, 0] + (vOffset - n * 2) * vScale, '--', 'LineWidth', dashedLineWidth, 'Color', 'r');
            plot([xWLocs(numQLeft+1), xWLocs(end)], [0, 0] + (vOffset - n * 2) * vScale, 'LineWidth', 2, 'Color', 'r');

            wXLocs = [xWLocs(1:numQLeft-1), xWLocs(end-numQRight+1:end)];
        end

        scatter(wXLocs, zeros(1, length(wXLocs)) + (vOffset - n * 2) * vScale,  400, 'r', '.');
        scatter(wXLocs(1) - 1, 0 + (vOffset - n * 2) * vScale, 80, 'r', 's', 'LineWidth', 1.5);
    end

    %% Draw texts
    textOffset = 0.2;

    if ~timeIndices
        idxP = "";
        idxV = "";
    else
        if n == 0
            idxP = "n+1";
            idxV = "n+1/2";
        else
            idxP = "n";
            idxV = "n-1/2";
        end
    end

    % pressures

    textPNext = [];
    textQNext = [];

    for idx = 0:numPLeft-1
        textPNext = [textPNext, "$p_{"+ num2str(idx) + "}^{"+ idxP + "}$"];
    end
    for idx = -numPRight+1:0
        if idx == 0
            num = "";
        else
            num = num2str(idx);
        end
        textPNext = [textPNext, "$p_{M_p" + num + "}^{"+ idxP + "}$"];
    end

    for idx = 0:numQLeft-1
        textQNext = [textQNext, "$q_{"+ num2str(idx) + "}^{"+ idxP + "}$"];
    end
    for idx = -numQRight+1:0
        if idx == 0
            num = "";
        else
            num = num2str(idx);
        end
        textQNext = [textQNext, "$q_{M_q" + num + "}^{"+ idxP + "}$"];
    end


    text(pXLocs, zeros(1, length(pXLocs)) + (textOffset - n * 2) * vScale, textPNext, 'interpreter', 'latex', ...
       'Fontsize', fontSize, 'color', 'b', 'horizontalAlignment', 'center');

    text(qXLocs, zeros(1, length(qXLocs)) + (textOffset - n * 2) * vScale, textQNext, 'interpreter', 'latex', ...
       'Fontsize', fontSize, 'color', 'r', 'horizontalAlignment', 'center');

    if n ~= 1 || drawVnmh
        % velocities

        textVNext = [];
        textWNext = [];

        if alignDashedToOutside
            rangeVEnd = numPLeft-2;                   
            rangeVStart = -numPRight+1;
            rangeWEnd = numQRight-1;
            rangeWStart = -numQLeft+2;
        else
            rangeVEnd = numPLeft-1;
            rangeVStart = -numPRight+2;
            rangeWEnd = numQRight-2;                   
            rangeWStart = -numQLeft+1;

        end
        for idx = 0:rangeVEnd
            textVNext = [textVNext, "$v_{"+ num2str(1 + idx*2) + "/2}^{"+ idxV + "}$"];
        end

        for idx = rangeVStart:0
            textVNext = [textVNext, "$v_{M_p" + num2str(-1 + idx*2) + "/2}^{"+ idxV + "}$"];
        end

        for idx = 0:rangeWEnd
            textWNext = [textWNext, "$w_{"+ num2str(1 + idx*2) + "/2}^{"+ idxV + "}$"];
        end
        for idx = rangeWStart:0
            textWNext = [textWNext, "$w_{M_q" + num2str(-1 + idx*2) + "/2}^{"+ idxV + "}$"];
        end

        text(vXLocs, zeros(1, length(vXLocs)) + (vOffset + textOffset - n * 2) * vScale, textVNext, 'interpreter', 'latex', ...
           'Fontsize', fontSize, 'color', 'b', 'horizontalAlignment', 'center', 'verticalAlignment', 'middle');


        text(wXLocs, zeros(1, length(wXLocs)) + (vOffset + textOffset - n * 2) * vScale, textWNext, 'interpreter', 'latex', ...
           'Fontsize', fontSize, 'color', 'r', 'horizontalAlignment', 'center', 'verticalAlignment', 'middle');

        vVelTextOffset = -0.8 * textOffset;

        text(vXLocs(end) + 1, (vOffset + vVelTextOffset - n * 2) * vScale, "$v_{M_p+1/2}^{"+ idxV + "}$", 'interpreter', 'latex', ...
           'Fontsize', fontSize, 'color', 'b', 'horizontalAlignment', 'center', 'verticalAlignment', 'middle');

        text(wXLocs(1) - 1, (vOffset + vVelTextOffset - n * 2) * vScale, "$w_{-1/2}^{"+ idxV + "}$", 'interpreter', 'latex', ...
           'Fontsize', fontSize, 'color', 'r', 'horizontalAlignment', 'center', 'verticalAlignment', 'middle');
    end

end
%% virtual points
vPressureTextOffset = -1.5 * textOffset;

pXMp1 = pXLocs(end) + 1;
scatter(pXMp1, 0 - n * 2 * vScale, 80, 'b', 'Linewidth', 2);
%             arrow([pXMp1, pXMp1], [vpTextOffset * 0.7, vpTextOffset * 0.2] - n * 2 * vScale, 1.5, 0.15, 0.25, 'b')

text(pXMp1, (vPressureTextOffset - n * 2) * vScale, "$p_{M_p+1}^{" + idxP + "}$", 'interpreter', 'latex', ...
   'Fontsize', fontSize, 'color', 'b', 'horizontalAlignment', 'center');

qXm1 = qXLocs(1) - 1;

scatter(qXm1, 0 - n * 2 * vScale, 80, 'r', 'Linewidth', 2);
text(qXm1, (vPressureTextOffset - n * 2) * vScale, "$q_{-1}^{" + idxP + "}$", 'interpreter', 'latex', ...
   'Fontsize', fontSize, 'color', 'r', 'horizontalAlignment', 'center');

%             arrow([qXm1, qXm1], [vpTextOffset * 0.7, vpTextOffset * 0.2] - n * 2 * vScale, 1.5, 0.15, 0.25, 'r')


set(gcf, 'color', 'w')
axis off


% plot options 
set(gca, 'Fontsize', fontSize)