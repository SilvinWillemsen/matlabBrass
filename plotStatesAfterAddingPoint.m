locsLeft = 0:length(up)-1;
locsRight = (0:length(wp)-1)+length(up) + alf;
if connectedWithP
    locsRight = locsRight - 1;
end
plotVel = false;


if ~plotVel
    %% Plot pressures

    subplot(2,1,1)
    if connectedWithP
        hold off;
        plot(locsLeft, up, '-o');
        hold on;
        plot(locsRight, wp, '-o');
    else
        hold off;
        plot([locsLeft, locsLeft(end) + 1], [up; upMp1], '-o');
        hold on;
        plot([locsRight(1) - 1, locsRight], [wpm1; wp], '-o');
        plot((locsLeft(end) + 1), upMp1, 'Marker', 'o', 'MarkerSize', 10, 'Color', 'r');
        plot((locsRight(1) - 1), wpm1, 'Marker', 'o', 'MarkerSize', 10,  'Color', 'b');

    end
    subplot(2,1,2)
    if connectedWithP
            hold off;
            plot(locsLeft, upPrev, '-o');
            hold on;
            plot(locsRight, wpPrev, '-o');
        else
            hold off;
            plot([locsLeft, locsLeft(end) + 1], [upPrev; upMp1Prev], '-o');
            hold on;
            plot([locsRight(1) - 1, locsRight], [wpm1; wpPrev], '-o');
            plot((locsLeft(end) + 1), upMp1Prev, 'Marker', 'o', 'MarkerSize', 10, 'Color', 'r');
            plot((locsRight(1) - 1), wpm1Prev, 'Marker', 'o', 'MarkerSize', 10,  'Color', 'b');
    end
else
    %         xlim([locsLeft(end-10), locsRight(10)])

    %% Plot velocities
    subplot(2,1,1)
    if connectedWithP
        hold off;
        plot(locsLeft + 0.5, [uvNext; uvNextMph], 'Marker', '.', 'MarkerSize', 10, 'Color', 'r');
        hold on;
        plot(locsRight - 0.5, [wvNextmh; wvNext], 'Marker', '.', 'MarkerSize', 10,  'Color', 'b');
        plot(locsLeft(end) + 0.5, uvNextMph, 'Marker', 'o', 'MarkerSize', 10, 'Color', 'r');
        plot(locsRight(1) - 0.5, wvNextmh, 'Marker', 'o', 'MarkerSize', 10,  'Color', 'b');
    else
        hold off;
        plot(locsLeft + 0.5, uvNext, 'Marker', '.', 'MarkerSize', 10, 'Color', 'r');
        hold on;
        plot(locsRight - 0.5, wvNext, 'Marker', '.', 'MarkerSize', 10,  'Color', 'b');

    end
    subplot(2,1,2)

 % plot prev states
    if connectedWithP
        hold off;
        plot(locsLeft + 0.5, [uv; uvMph], 'Marker', '.', 'MarkerSize', 10, 'Color', 'r');
        hold on;
        plot(locsRight - 0.5, [wvmh; wv], 'Marker', '.', 'MarkerSize', 10,  'Color', 'b');
        plot(locsLeft(end) + 0.5, uvMph, 'Marker', 'o', 'MarkerSize', 10, 'Color', 'r');
        plot(locsRight(1) - 0.5, wvmh, 'Marker', 'o', 'MarkerSize', 10,  'Color', 'b');
    else
        hold off;
        plot(locsLeft + 0.5, uv, 'Marker', '.', 'MarkerSize', 10, 'Color', 'r');
        hold on;
        plot(locsRight - 0.5, wv, 'Marker', '.', 'MarkerSize', 10,  'Color', 'b');

    end
end
%         xlim([locsLeft(end-10), locsRight(10)])
pause(0.1)
drawnow;